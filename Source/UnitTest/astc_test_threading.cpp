// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2026 Arm Limited
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not
// use this file except in compliance with the License. You may obtain a copy
// of the License at:
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
// WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
// License for the specific language governing permissions and limitations
// under the License.
// ----------------------------------------------------------------------------

/**
 * @file
 * @brief Randomized public API threading and context lifetime test.
 *
 * Generate a fixed corpus of 2D and 3D images and record single-threaded compressed and
 * decompressed reference bytes for each compatible image/settings pair. Settings cover multiple
 * profiles, block sizes, quality presets, and flags. Every completed compression and decompression
 * in the main test loop must match its reference byte-for-byte.
 *
 * Each worker runs an independent randomized loop, initializing configurations or creating and
 * reusing its own active context. Configuration initialization mixes profiles, block sizes,
 * quality presets, and flags, and must match the serial configuration field-for-field.
 * Contexts are configured to use between one and --max-threads codec threads; each operation uses
 * a randomized subset of distinct thread indices with staggered starts. --workers controls the
 * number of independent loops, so total codec concurrency can reach workers * max-threads.
 *
 * Bounded shared pools let workers create children from other workers' contexts, including active
 * parents, and decompress completed images compressed by other workers. Context chains retain
 * their root while allowing intermediate parents to be freed. A context processes only one image
 * at a time, and all participating codec threads join before the context is reset or freed.
 *
 * Compression cancellation is requested externally or from progress callbacks. Cancelled output is
 * discarded, then the same context is reset and used for a complete compression checked against
 * its reference. Progress callbacks are also checked for serialization, monotonicity, and a final
 * 100% report on complete compressions.
 *
 * Workers stop issuing jobs after --duration seconds. A watchdog checks total runtime, including
 * reference generation and worker shutdown, and fails the process after the requested duration
 * plus 30 seconds.
 *
 * Aggregate operation counts are logged every five seconds and at shutdown. Final counters report
 * coverage for each image/settings pair; failures identify the seed, worker, iteration, and
 * reference. Seeds control random choices, but scheduling and shared-pool contents make exact
 * replay nondeterministic. Run under ASAN or TSAN to detect memory errors and data races that
 * output comparisons may not reveal.
 */

#include "astcenc.h"

#include <algorithm>
#include <atomic>
#include <barrier>
#include <charconv>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <span>
#include <string_view>
#include <thread>
#include <tuple>
#include <vector>

using Clock = std::chrono::steady_clock;

/** @brief Command line settings. */
struct Options
{
	/** @brief Duration of the randomized loop, in seconds. */
	unsigned int duration { 30 };

	/** @brief Number of independent workers. */
	unsigned int workers { 4 };

	/** @brief Maximum number of threads configured per context. */
	unsigned int max_threads { 4 };

	/** @brief Random seed; scheduling remains nondeterministic. */
	unsigned int seed { 1 };
};

/** @brief Fixed image data, available in both LDR and HDR representations. */
struct Image
{
	/** @brief Image width. */
	size_t x;

	/** @brief Image height. */
	size_t y;

	/** @brief Image depth. */
	size_t z;

	/** @brief LDR RGBA texels. */
	std::vector<uint8_t> ldr;

	/** @brief HDR RGBA texels. */
	std::vector<float> hdr;

	/**
	 * @brief Create a view whose plane array is owned by the caller.
	 *
	 * The image storage and plane array must remain valid while the view is in use.
	 *
	 * @param[in]  is_hdr   True to use HDR float texels, false to use LDR byte texels.
	 * @param[out] planes   Plane pointers populated for the selected representation.
	 *
	 * @return A view borrowing this image's texel storage and the caller's plane array.
	 */
	astcenc_image view(bool is_hdr, std::vector<void*>& planes)
	{
		planes.clear();
		for (size_t i = 0; i < z; i++)
		{
			if (is_hdr)
			{
				planes.push_back(hdr.data() + i * x * y * 4);
			}
			else
			{
				planes.push_back(ldr.data() + i * x * y * 4);
			}
		}

		astcenc_type type = ASTCENC_TYPE_U8;
		if (is_hdr)
		{
			type = ASTCENC_TYPE_F32;
		}

		return astcenc_image {
			.dim_x = x,
			.dim_y = y,
			.dim_z = z,
			.data_type = type,
			.data = planes.data()
		};
	}

	/**
	 * @brief View the selected representation as bytes for exact comparison.
	 *
	 * @param is_hdr   True to view HDR float texels, false to view LDR byte texels.
	 *
	 * @return A read-only byte view valid until the selected storage is freed or reallocated.
	 */
	std::span<const std::byte> as_bytes(bool is_hdr) const
	{
		if (is_hdr)
		{
			return std::as_bytes(std::span(hdr));
		}

		return std::as_bytes(std::span(ldr));
	}
};

/** @brief One configuration and its compatible image references. */
struct Setting
{
	/** @brief Quality preset used to initialize the configuration. */
	float quality;

	/** @brief Unmodified configuration initialized serially for comparison. */
	astcenc_config initialized_config;

	/** @brief Configuration inherited by child contexts. */
	astcenc_config config;

	/** @brief True if input and output use floating point texels. */
	bool is_hdr;

	/** @brief Indices of compatible reference images. */
	std::vector<size_t> references;
};

/** @brief Expected bytes for one image and configuration. */
struct Reference
{
	/** @brief Configuration index. */
	size_t setting;

	/** @brief Fixed image index. */
	size_t image;

	/** @brief Expected compressed bytes. */
	std::vector<uint8_t> compressed;

	/** @brief Expected decompressed bytes. */
	std::vector<std::byte> decompressed;
};

/** @brief Context ownership that keeps the root alive for all descendants. */
struct Context
{
	/** @brief Codec context. */
	astcenc_context* codec { nullptr };

	/** @brief Root owner; intermediate parents can be freed independently. */
	std::shared_ptr<Context> root;

	/** @brief Configuration index. */
	size_t setting;

	/** @brief Configured thread count. */
	size_t threads;

	/** @brief Distance from the root context. */
	size_t depth { 0 };

	/** @brief Free this context before releasing its root owner. */
	~Context()
	{
		astcenc_context_free(codec);
	}
};

/** @brief Completed compression passed to a different worker. */
struct Packet
{
	/** @brief Reference index. */
	size_t reference;

	/** @brief Worker that produced these bytes. */
	size_t producer;

	/** @brief Actual compressed bytes from the producer. */
	std::vector<uint8_t> data;
};

/** @brief Atomic diagnostic snapshot for a worker that may become stuck. */
struct WorkerStatus
{
	/** @brief Most recent job iteration. */
	std::atomic<size_t> iteration { 0 };

	/** @brief Most recent reference index. */
	std::atomic<size_t> reference { 0 };

	/** @brief Current operation; points to a string literal. */
	std::atomic<const char*> operation { "Starting" };
};

/** @brief Shared corpus, bounded exchange pools, and coverage counters. */
struct State
{
	/** @brief Run settings. */
	Options options;

	/** @brief Immutable input corpus after reference generation. */
	std::vector<Image> images;

	/** @brief Immutable configurations after reference generation. */
	std::vector<Setting> settings;

	/** @brief Immutable serial results after reference generation. */
	std::vector<Reference> references;

	/** @brief Per-reference compression counts, protected by @c lock. */
	std::vector<size_t> compression_hits;

	/** @brief Per-reference decompression counts, protected by @c lock. */
	std::vector<size_t> decompression_hits;

	/** @brief Per-worker diagnostic snapshots for the watchdog. */
	std::unique_ptr<WorkerStatus[]> workers;

	/** @brief Protects the exchange pools and failure diagnostics. */
	std::mutex lock;

	/** @brief Contexts available as parents, including contexts currently in use. */
	std::vector<std::shared_ptr<Context>> parents;

	/** @brief Completed images available for another worker to decode. */
	std::vector<Packet> packets;

	/** @brief True after the first failure. */
	std::atomic<bool> is_failed { false };

	/** @brief Number of configuration initializations compared with their serial references. */
	std::atomic<size_t> configurations { 0 };

	/** @brief Number of standalone contexts created. */
	std::atomic<size_t> roots { 0 };

	/** @brief Number of child contexts created. */
	std::atomic<size_t> children { 0 };

	/** @brief Number of descendants created from non-root parents. */
	std::atomic<size_t> grandchildren { 0 };

	/** @brief Number of single-threaded contexts created. */
	std::atomic<size_t> single { 0 };

	/** @brief Number of multi-threaded contexts created. */
	std::atomic<size_t> multi { 0 };

	/** @brief Number of completed compressions compared with their references. */
	std::atomic<size_t> compressions { 0 };

	/** @brief Number of completed decompressions compared with their references. */
	std::atomic<size_t> decompressions { 0 };

	/** @brief Number of compression cancellation requests issued externally. */
	std::atomic<size_t> cancellations { 0 };

	/** @brief Number of compression cancellations issued by callbacks. */
	std::atomic<size_t> callback_cancellations { 0 };

	/** @brief Number of images decoded by a different worker than their producer. */
	std::atomic<size_t> transfers { 0 };
};

/**
 * @brief Print a cumulative snapshot of aggregate operation counts.
 *
 * @param state     The test state.
 * @param elapsed   The elapsed run time, in seconds.
 */
static void print_statistics(const State& state, double elapsed)
{
	std::printf("\nStatistics (%.1f s):\n"
	            "  - Configurations validated:       %8zu\n"
	            "  - Root contexts:                  %8zu\n"
	            "  - Child contexts:                 %8zu\n"
	            "  - Grandchild contexts:            %8zu\n"
	            "  - Single-threaded contexts:       %8zu\n"
	            "  - Multi-threaded contexts:        %8zu\n"
	            "  - Compressions validated:         %8zu\n"
	            "  - Decompressions validated:       %8zu\n"
	            "  - Cross-worker transfers:         %8zu\n"
	            "  - External cancellation requests: %8zu\n"
	            "  - Callback cancellation requests: %8zu\n",
	            elapsed, state.configurations.load(),
	            state.roots.load(), state.children.load(), state.grandchildren.load(),
	            state.single.load(), state.multi.load(), state.compressions.load(),
	            state.decompressions.load(), state.transfers.load(),
	            state.cancellations.load(), state.callback_cancellations.load());

	std::fflush(stdout);
}

/**
 * @brief Count the decimal digits needed to print an unsigned value.
 *
 * @param value   The value to measure.
 *
 * @return The required column width, including one digit for zero.
 */
static int decimal_width(size_t value)
{
	int width = 1;
	while (value >= 10)
	{
		value /= 10;
		width++;
	}

	return width;
}

/**
 * @brief Print per-reference operation counts with columns sized to their contents.
 *
 * @param state   The test state after all workers have stopped.
 */
static void print_references(const State& state)
{
	int reference_width = 1;
	int setting_width = 1;
	int image_width = 1;
	int comp_width = 1;
	int decomp_width = 1;

	for (size_t i = 0; i < state.references.size(); i++)
	{
		const Reference& reference = state.references[i];
		reference_width = std::max(reference_width, decimal_width(i));
		setting_width = std::max(setting_width, decimal_width(reference.setting));
		image_width = std::max(image_width, decimal_width(reference.image));
		comp_width = std::max(comp_width, decimal_width(state.compression_hits[i]));
		decomp_width = std::max(decomp_width, decimal_width(state.decompression_hits[i]));
	}

	std::puts("\nReferences:");
	for (size_t i = 0; i < state.references.size(); i++)
	{
		const Reference& reference = state.references[i];
		std::printf("  - %*zu (opt=%*zu img=%*zu):   comp=%*zu   decomp=%*zu\n",
		            reference_width, i,
		            setting_width, reference.setting,
		            image_width, reference.image,
		            comp_width, state.compression_hits[i],
		            decomp_width, state.decompression_hits[i]);
	}
}

/** @brief Per-operation progress state, independent of the library locks. */
struct Progress
{
	/** @brief Context being compressed. */
	astcenc_context* context;

	/** @brief True if the callback should request cancellation. */
	bool is_cancel_enabled;

	/** @brief Minimum progress at which to request cancellation. */
	float cancel_at;

	/** @brief Detect overlapping callbacks without introducing a test data race. */
	std::atomic_flag is_active = ATOMIC_FLAG_INIT;

	/** @brief Last reported progress. */
	std::atomic<float> last { 0.0f };

	/** @brief True if progress was invalid or callbacks overlapped. */
	std::atomic<bool> is_invalid { false };

	/** @brief True if cancellation was requested by a callback. */
	std::atomic<bool> is_cancelled { false };
};

/** @brief Operation state for callbacks invoked on the calling codec thread. */
static thread_local Progress* current_progress;

/**
 * @brief Validate progress and optionally cancel from inside the callback.
 *
 * @param value   The reported progress percentage.
 */
static void progress_callback(float value)
{
	Progress& progress = *current_progress;
	if (progress.is_active.test_and_set())
	{
		progress.is_invalid = true;
		return;
	}

	float last = progress.last.exchange(value);
	if (!std::isfinite(value) || (value < last) || (value < 0.0f) || (value > 100.0f))
	{
		progress.is_invalid = true;
	}

	if (progress.is_cancel_enabled && (value >= progress.cancel_at) && !progress.is_cancelled.exchange(true))
	{
		astcenc_error status = astcenc_compress_cancel(progress.context);
		if (status != ASTCENC_SUCCESS)
		{
			progress.is_invalid = true;
		}
	}

	progress.is_active.clear();
}

/**
 * @brief Report enough information to identify the failing randomized job.
 *
 * @param[in,out] state       Shared state to mark as failed.
 * @param         worker      Worker index for diagnostics.
 * @param         iteration   Job iteration for diagnostics.
 * @param         reference   Reference index for diagnostics.
 * @param[in]     message     Failure description to log.
 *
 * @return Always @c false, allowing callers to propagate failure directly.
 */
static bool fail(
	State& state,
	size_t worker,
	size_t iteration,
	size_t reference,
	const char* message)
{
	std::lock_guard<std::mutex> guard(state.lock);

	state.is_failed = true;
	std::fprintf(stderr, "\nFAIL seed=%u worker=%zu iteration=%zu reference=%zu: %s\n",
	             state.options.seed, worker, iteration, reference, message);
	return false;
}

/**
 * @brief Select configuration fields for comparison without including structure padding.
 *
 * @param config   The configuration whose fields are selected.
 *
 * @return A tuple of field references, valid while @p config remains valid.
 */
static auto config_fields(const astcenc_config& config)
{
	return std::tie(
		config.profile,
		config.flags,
		config.block_x,
		config.block_y,
		config.block_z,
		config.cw_r_weight,
		config.cw_g_weight,
		config.cw_b_weight,
		config.cw_a_weight,
		config.a_scale_radius,
		config.rgbm_m_scale,
		config.tune_partition_count_limit,
		config.tune_2partition_index_limit,
		config.tune_3partition_index_limit,
		config.tune_4partition_index_limit,
		config.tune_block_mode_limit,
		config.tune_refinement_limit,
		config.tune_candidate_limit,
		config.tune_2partitioning_candidate_limit,
		config.tune_3partitioning_candidate_limit,
		config.tune_4partitioning_candidate_limit,
		config.tune_db_limit,
		config.tune_mse_overshoot,
		config.tune_2partition_early_out_limit_factor,
		config.tune_3partition_early_out_limit_factor,
		config.tune_2plane_early_out_limit_correlation,
		config.tune_search_mode0_enable,
		config.progress_callback);
}

/**
 * @brief Initialize a configuration and compare it with its serial reference.
 *
 * @param[in,out] state       Shared configurations, counters, and failure state.
 * @param         setting     Configuration index to initialize.
 * @param         worker      Worker index for diagnostics.
 * @param         iteration   Job iteration for diagnostics.
 *
 * @return @c true if initialization matches the reference, @c false on failure.
 */
static bool initialize_config(State& state, size_t setting, size_t worker, size_t iteration)
{
	const Setting& expected = state.settings[setting];
	const astcenc_config& reference = expected.initialized_config;
	size_t reference_index = expected.references.front();
	state.workers[worker].reference = reference_index;
	state.workers[worker].operation = "Initializing config";

	astcenc_config config;
	astcenc_error status = astcenc_config_init(reference.profile, reference.block_x, reference.block_y,
	                                           reference.block_z, expected.quality, reference.flags, &config);
	if (status != ASTCENC_SUCCESS)
	{
		return fail(state, worker, iteration, reference_index, astcenc_get_error_string(status));
	}

	if (config_fields(config) != config_fields(reference))
	{
		return fail(state, worker, iteration, reference_index, "Config fields differ from serial reference");
	}

	state.configurations++;
	return true;
}

/**
 * @brief Create a standalone context or a child of a shared, possibly active parent.
 *
 * @param[in,out] state       Shared configurations, counters, and failure state.
 * @param         setting     Configuration index for the new context.
 * @param         threads     Number of codec threads to configure.
 * @param[in]     parent      Optional parent context with the same configuration.
 * @param         worker      Worker index for diagnostics.
 * @param         iteration   Job iteration for diagnostics.
 * @param         reference   Reference index for diagnostics.
 *
 * @return The new context, or @c nullptr if codec context allocation fails.
 */
static std::shared_ptr<Context> make_context(
	State& state,
	size_t setting,
	size_t threads,
	const std::shared_ptr<Context>& parent,
	size_t worker,
	size_t iteration,
	size_t reference)
{
	auto context = std::make_shared<Context>();
	context->setting = setting;
	context->threads = threads;
	astcenc_error status;
	if (parent)
	{
		// Keep the root alive without retaining intermediate parents
		context->root = parent;
		if (parent->root)
		{
			context->root = parent->root;
		}

		context->depth = parent->depth + 1;
		status = astcenc_context_alloc(nullptr, threads, &context->codec, parent->codec);
	}
	else
	{
		status = astcenc_context_alloc(&state.settings[setting].config, threads, &context->codec, nullptr);
	}

	if (status != ASTCENC_SUCCESS)
	{
		fail(state, worker, iteration, reference, astcenc_get_error_string(status));
		return nullptr;
	}

	if (parent)
	{
		state.children++;
		if (parent->depth > 0)
		{
			state.grandchildren++;
		}
	}
	else
	{
		state.roots++;
	}

	if (threads == 1)
	{
		state.single++;
	}
	else
	{
		state.multi++;
	}

	return context;
}

/**
 * @brief Populate a small fixed corpus with smooth, noisy, and transparent regions.
 *
 * @param[in,out] state   Shared state to populate with fixed input images.
 */
static void make_images(State& state)
{
	const size_t dimensions[][3] { { 33, 29, 1 }, { 127, 65, 1 }, { 257, 513, 1 }, { 25, 19, 5 } };
	std::mt19937 random(0xA57C);
	for (const auto& dim : dimensions)
	{
		Image image { dim[0], dim[1], dim[2], {}, {} };
		size_t count = image.x * image.y * image.z * 4;
		image.ldr.resize(count);
		image.hdr.resize(count);

		// Mix gradients, noise, and zero alpha in the same reusable image
		for (size_t i = 0; i < count; i++)
		{
			unsigned int value = static_cast<unsigned int>((i / 4) % 256);
			if (((i / (image.x * 4)) % 3) == 0)
			{
				value = random() % 256;
			}

			if (((i % 4) == 3) && (((i / 4) % 7) == 0))
			{
				value = 0;
			}

			image.ldr[i] = static_cast<uint8_t>(value);
			image.hdr[i] = static_cast<float>(value) / 255.0f;
			if ((i % 4) != 3)
			{
				image.hdr[i] *= 8.0f;
			}
		}

		state.images.push_back(std::move(image));
	}
}

/**
 * @brief Build serial byte references for every compatible image and settings pair.
 *
 * @param[in,out] state   Shared corpus, configurations, and reference results to populate.
 *
 * @return @c true on success, @c false if a codec API call fails.
 */
static bool make_references(State& state)
{
	// Cover different profiles, block sizes, quality levels, and preprocessing paths
	const astcenc_profile profiles[] {
		ASTCENC_PRF_LDR, ASTCENC_PRF_LDR, ASTCENC_PRF_LDR_SRGB, ASTCENC_PRF_LDR,
		ASTCENC_PRF_HDR, ASTCENC_PRF_HDR_RGB_LDR_A, ASTCENC_PRF_LDR, ASTCENC_PRF_LDR
	};
	const unsigned int blocks[] { 4, 6, 8, 4, 6, 4, 4, 6 };
	const float qualities[] {
		ASTCENC_PRE_FASTEST, ASTCENC_PRE_FAST, ASTCENC_PRE_MEDIUM, ASTCENC_PRE_MEDIUM,
		ASTCENC_PRE_FAST, ASTCENC_PRE_FASTEST, ASTCENC_PRE_FASTEST, ASTCENC_PRE_FAST
	};
	const unsigned int flags[] {
		0, ASTCENC_FLG_USE_ALPHA_WEIGHT, ASTCENC_FLG_USE_DECODE_UNORM8, ASTCENC_FLG_MAP_NORMAL,
		0, 0, 0, ASTCENC_FLG_DECOMPRESS_ONLY
	};
	const astcenc_swizzle swizzle { ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A };

	for (size_t s = 0; s < std::size(profiles); s++)
	{
		unsigned int block_z = 1;
		if (s == 6)
		{
			block_z = 4;
		}

		Setting setting {};
		setting.quality = qualities[s];
		setting.is_hdr = (profiles[s] == ASTCENC_PRF_HDR) || (profiles[s] == ASTCENC_PRF_HDR_RGB_LDR_A);
		astcenc_error status = astcenc_config_init(profiles[s], blocks[s], blocks[s], block_z,
		                                           setting.quality, flags[s], &setting.initialized_config);
		if (status != ASTCENC_SUCCESS)
		{
			return fail(state, 0, 0, s, astcenc_get_error_string(status));
		}

		setting.config = setting.initialized_config;
		setting.config.progress_callback = progress_callback;
		if (s == 1)
		{
			setting.config.a_scale_radius = 4;
		}

		state.settings.push_back(std::move(setting));

		// Reference generation also needs compression for decompress-only test contexts
		astcenc_config reference_config = state.settings[s].config;
		reference_config.flags &= ~ASTCENC_FLG_DECOMPRESS_ONLY;
		reference_config.progress_callback = nullptr;
		astcenc_context* codec = nullptr;
		status = astcenc_context_alloc(&reference_config, 1, &codec, nullptr);
		if (status != ASTCENC_SUCCESS)
		{
			return fail(state, 0, 0, s, astcenc_get_error_string(status));
		}

		std::unique_ptr<astcenc_context, decltype(&astcenc_context_free)> owner(codec, astcenc_context_free);
		for (size_t i = 0; i < state.images.size(); i++)
		{
			Image image = state.images[i];
			if ((image.z > 1) != (block_z > 1))
			{
				continue;
			}

			// Record both stages serially before any workers share contexts or compressed data
			std::vector<void*> planes;
			astcenc_image view = image.view(state.settings[s].is_hdr, planes);
			size_t block_count = ((image.x + blocks[s] - 1) / blocks[s])
			                   * ((image.y + blocks[s] - 1) / blocks[s])
			                   * ((image.z + block_z - 1) / block_z);
			Reference reference { s, i, std::vector<uint8_t>(block_count * 16), {} };
			status = astcenc_compress_image(codec, &view, &swizzle,
			                                 reference.compressed.data(), reference.compressed.size(), 0);
			if (status == ASTCENC_SUCCESS)
			{
				status = astcenc_decompress_image(codec, reference.compressed.data(),
				                                   reference.compressed.size(), &view, &swizzle, 0);
			}

			if (status != ASTCENC_SUCCESS)
			{
				return fail(state, 0, 0, s, astcenc_get_error_string(status));
			}

			// Own the reference bytes because this temporary image is about to be destroyed
			auto bytes = image.as_bytes(state.settings[s].is_hdr);
			reference.decompressed.assign(bytes.begin(), bytes.end());
			state.settings[s].references.push_back(state.references.size());
			state.references.push_back(std::move(reference));
			state.compression_hits.push_back(0);
			state.decompression_hits.push_back(0);
		}
	}

	return true;
}

/** @brief How a compression should be cancelled. */
enum class Cancel
{
	NONE,
	EXTERNAL,
	CALLBACK
};

/**
 * @brief Run an operation using a randomized subset of distinct codec thread indices.
 *
 * All participating codec threads join before returning. Cancelled compression output is undefined
 * and must not be compared with its reference.
 *
 * @param[in,out] state            Shared counters and failure diagnostics.
 * @param[in,out] context          Context reserved for this operation.
 * @param[in,out] image            Input image for compression, output image for decompression.
 * @param[in,out] data             Output blocks for compression, input blocks for decompression.
 * @param         is_compression   True to compress, false to decompress.
 * @param         cancel           Cancellation mode; must be @c Cancel::NONE for decompression.
 * @param[in,out] random           Generator for codec thread selection and scheduling.
 * @param         worker           Worker index for diagnostics.
 * @param         iteration        Job iteration for diagnostics.
 * @param         reference        Reference index for diagnostics.
 *
 * @return @c true if codec calls and progress checks succeed, @c false on failure.
 */
static bool run_codec(
	State& state,
	Context& context,
	Image& image,
	std::vector<uint8_t>& data,
	bool is_compression,
	Cancel cancel,
	std::mt19937& random,
	size_t worker,
	size_t iteration,
	size_t reference)
{
	const char* operation = "Decompressing";
	if (is_compression)
	{
		operation = "Compressing";
		if (cancel != Cancel::NONE)
		{
			operation = "Compressing with cancellation";
		}
	}

	state.workers[worker].operation = operation;

	// Reusing a context must reset all state left by its previous operation
	astcenc_error status;
	if (is_compression)
	{
		status = astcenc_compress_reset(context.codec);
	}
	else
	{
		status = astcenc_decompress_reset(context.codec);
	}

	if (status != ASTCENC_SUCCESS)
	{
		return fail(state, worker, iteration, reference, astcenc_get_error_string(status));
	}

	std::vector<void*> planes;
	astcenc_image view = image.view(state.settings[context.setting].is_hdr, planes);
	const astcenc_swizzle swizzle { ASTCENC_SWZ_R, ASTCENC_SWZ_G, ASTCENC_SWZ_B, ASTCENC_SWZ_A };
	Progress progress {
		.context = context.codec,
		.is_cancel_enabled = cancel == Cancel::CALLBACK,
		.cancel_at = 1.0f + static_cast<float>(random() % 100)
	};

	// Use distinct thread indices and vary how many configured threads actually participate
	std::vector<size_t> indices(context.threads);
	std::iota(indices.begin(), indices.end(), 0);
	std::shuffle(indices.begin(), indices.end(), random);
	size_t active = 1 + random() % context.threads;
	std::vector<astcenc_error> results(active);
	std::vector<std::thread> threads;
	std::barrier start(static_cast<std::ptrdiff_t>(active + 1));

	for (size_t i = 0; i < active; i++)
	{
		unsigned int delay = random() % 100;
		if ((cancel == Cancel::EXTERNAL) && (i > 0))
		{
			// Exercise cancellation while some participating threads have not entered the codec
			delay += 1000;
		}

		threads.emplace_back([&, i, delay] {
			current_progress = &progress;
			start.arrive_and_wait();
			std::this_thread::sleep_for(std::chrono::microseconds(delay));
			if (is_compression)
			{
				results[i] = astcenc_compress_image(context.codec, &view, &swizzle,
				                                     data.data(), data.size(), indices[i]);
			}
			else
			{
				results[i] = astcenc_decompress_image(context.codec, data.data(), data.size(),
				                                       &view, &swizzle, indices[i]);
			}

			current_progress = nullptr;
		});
	}

	start.arrive_and_wait();
	if (cancel == Cancel::EXTERNAL)
	{
		std::this_thread::sleep_for(std::chrono::microseconds(random() % 500));
		status = astcenc_compress_cancel(context.codec);
		state.cancellations++;
	}

	// Drain every codec call before checking results or allowing the context to be reused
	for (auto& thread : threads)
	{
		thread.join();
	}

	if (status != ASTCENC_SUCCESS)
	{
		return fail(state, worker, iteration, reference, astcenc_get_error_string(status));
	}

	for (astcenc_error result : results)
	{
		if (result != ASTCENC_SUCCESS)
		{
			return fail(state, worker, iteration, reference, astcenc_get_error_string(result));
		}
	}

	if (progress.is_cancelled)
	{
		state.callback_cancellations++;
	}

	if (progress.is_invalid || (is_compression && (cancel == Cancel::NONE) && (progress.last != 100.0f)))
	{
		return fail(state, worker, iteration, reference, "Invalid or incomplete progress callbacks");
	}

	return true;
}

/**
 * @brief Mix configuration initialization, context lifetimes, and codec operations.
 *
 * @param[in,out] state      Shared corpus, exchange pools, counters, and failure state.
 * @param         worker     Index of this independent worker.
 * @param         deadline   Time after which no new jobs should start.
 */
static void run_worker(State& state, size_t worker, Clock::time_point deadline)
{
	std::seed_seq seeds { state.options.seed, static_cast<unsigned int>(worker) };
	std::mt19937 random(seeds);
	std::shared_ptr<Context> context;
	for (size_t iteration = 0; (Clock::now() < deadline) && !state.is_failed; iteration++)
	{
		state.workers[worker].iteration = iteration;
		size_t setting = random() % state.settings.size();
		if ((random() % 4) == 0)
		{
			bool is_successful = initialize_config(state, setting, worker, iteration);
			if (!is_successful)
			{
				return;
			}

			continue;
		}

		// Retain some contexts across jobs while still varying settings and reference images
		if (context && ((random() % 2) == 0))
		{
			setting = context->setting;
		}

		const auto& choices = state.settings[setting].references;
		size_t reference_index = choices[random() % choices.size()];
		bool is_compression = (random() % 3) != 0;
		if (state.settings[setting].config.flags & ASTCENC_FLG_DECOMPRESS_ONLY)
		{
			is_compression = false;
		}

		// Prefer actual compressed output from a different worker when decoding
		Packet packet {};
		bool is_transferred = false;
		if (!is_compression && !(state.settings[setting].config.flags & ASTCENC_FLG_DECOMPRESS_ONLY))
		{
			std::lock_guard<std::mutex> guard(state.lock);

			for (size_t i = 0; i < state.packets.size(); i++)
			{
				if (state.packets[i].producer != worker)
				{
					packet = std::move(state.packets[i]);
					state.packets.erase(state.packets.begin() + static_cast<std::ptrdiff_t>(i));
					reference_index = packet.reference;
					setting = state.references[reference_index].setting;
					is_transferred = true;
					break;
				}
			}
		}

		state.workers[worker].reference = reference_index;
		if (!context || (context->setting != setting) || ((random() % 3) == 0))
		{
			state.workers[worker].operation = "Allocating context";
			std::shared_ptr<Context> parent;
			std::shared_ptr<Context> evicted;
			{
				std::lock_guard<std::mutex> guard(state.lock);

				if (!state.parents.empty() && ((random() % 3) != 0))
				{
					size_t offset = random() % state.parents.size();
					for (size_t i = 0; i < state.parents.size(); i++)
					{
						auto& candidate = state.parents[(offset + i) % state.parents.size()];
						if (candidate->setting == setting)
						{
							parent = candidate;
							break;
						}
					}
				}
			}

			// Release the pool lock before allocation; the selected parent may still be active
			size_t threads = 1 + random() % state.options.max_threads;
			context = make_context(state, setting, threads, parent, worker, iteration, reference_index);
			if (!context)
			{
				return;
			}

			{
				std::lock_guard<std::mutex> guard(state.lock);

				// Bound retained contexts and release evicted owners after dropping the pool lock
				if (state.parents.size() < 24)
				{
					state.parents.push_back(context);
				}
				else
				{
					size_t slot = random() % state.parents.size();
					evicted = std::move(state.parents[slot]);
					state.parents[slot] = context;
				}
			}
		}

		const Reference& reference = state.references[reference_index];
		Image image = state.images[reference.image];
		std::vector<uint8_t> data;
		if (is_compression)
		{
			data.resize(reference.compressed.size(), 0xCD);
			Cancel cancel = static_cast<Cancel>(random() % 3);
			if (cancel != Cancel::NONE)
			{
				bool is_successful = run_codec(state, *context, image, data, true, cancel,
				                               random, worker, iteration, reference_index);
				if (!is_successful)
				{
					return;
				}

				// Cancelled data is undefined; verify recovery on the same context after reset
				std::fill(data.begin(), data.end(), 0xCD);
			}

			bool is_successful = run_codec(state, *context, image, data, true, Cancel::NONE,
			                               random, worker, iteration, reference_index);
			if (!is_successful)
			{
				return;
			}

			if (data != reference.compressed)
			{
				fail(state, worker, iteration, reference_index, "Compressed bytes differ from serial reference");
				return;
			}

			state.compressions++;
			std::lock_guard<std::mutex> guard(state.lock);

			state.compression_hits[reference_index]++;
			if (state.packets.size() == 48)
			{
				state.packets.erase(state.packets.begin());
			}

			state.packets.push_back(Packet { reference_index, worker, std::move(data) });
		}
		else
		{
			data = std::move(packet.data);
			if (!is_transferred)
			{
				data = reference.compressed;
			}

			// Poison the output so incomplete writes cannot match the serial reference
			std::fill(image.ldr.begin(), image.ldr.end(), 0xCD);
			std::fill(image.hdr.begin(), image.hdr.end(), -1234.0f);
			bool is_successful = run_codec(state, *context, image, data, false, Cancel::NONE,
			                               random, worker, iteration, reference_index);
			if (!is_successful)
			{
				return;
			}

			auto bytes = image.as_bytes(state.settings[setting].is_hdr);
			if (!std::ranges::equal(bytes, reference.decompressed))
			{
				fail(state, worker, iteration, reference_index, "Decompressed bytes differ from serial reference");
				return;
			}

			state.decompressions++;
			{
				std::lock_guard<std::mutex> guard(state.lock);

				state.decompression_hits[reference_index]++;
			}
			if (is_transferred)
			{
				state.transfers++;
			}
		}
	}

	state.workers[worker].operation = "Finished";
}

/**
 * @brief Parse bounded numeric command line options.
 *
 * @param         argc      Number of command line arguments, including the executable name.
 * @param[in]     argv      Command line argument strings.
 * @param[in,out] options   Defaults to update with parsed argument values.
 *
 * @return @c true for valid arguments, @c false for help or invalid arguments.
 */
static bool parse_options(int argc, char** argv, Options& options)
{
	for (int i = 1; i < argc; i += 2)
	{
		std::string_view name(argv[i]);
		if (name == "--help")
		{
			return false;
		}

		unsigned int* target = nullptr;
		unsigned int minimum = 0;
		unsigned int maximum = std::numeric_limits<unsigned int>::max();
		if (name == "--duration")
		{
			target = &options.duration;
			minimum = 1;
			maximum = 86400;
		}
		else if (name == "--workers")
		{
			target = &options.workers;
			minimum = 2;
			maximum = 64;
		}
		else if (name == "--max-threads")
		{
			target = &options.max_threads;
			minimum = 1;
			maximum = 16;
		}
		else if (name == "--seed")
		{
			target = &options.seed;
		}
		else
		{
			std::fprintf(stderr, "ERROR: Unknown option '%s'\n", argv[i]);
			return false;
		}

		if (i + 1 == argc)
		{
			std::fprintf(stderr, "ERROR: Missing value for %s\n", argv[i]);
			return false;
		}

		std::string_view text(argv[i + 1]);
		unsigned int value = 0;
		auto result = std::from_chars(text.data(), text.data() + text.size(), value);
		if ((result.ec != std::errc()) || (result.ptr != text.data() + text.size()) ||
		    (value < minimum) || (value > maximum))
		{
			std::fprintf(stderr, "ERROR: %s requires an unsigned integer from %u to %u (got '%s')\n",
			             argv[i], minimum, maximum, argv[i + 1]);
			return false;
		}

		*target = value;
	}

	return true;
}

/**
 * @brief Run the test with a deadline and a watchdog for deadlocked operations.
 *
 * @param     argc   Number of command line arguments, including the executable name.
 * @param[in] argv   Command line argument strings.
 *
 * @return 0 if the test passes or help is requested, 1 on failure.
 */
int main(int argc, char** argv)
{
	State state {};
	if (!parse_options(argc, argv, state.options))
	{
		std::fprintf(stderr, "Usage: %s [--duration seconds] [--workers count] "
		             "[--max-threads count] [--seed value]\n"
		             "Defaults: --duration 30 --workers 4 --max-threads 4 --seed 1\n"
		             "Ranges: --duration 1..86400, --workers 2..64, "
		             "--max-threads 1..16, --seed 0..%u\n",
		             argv[0], std::numeric_limits<unsigned int>::max());
		if ((argc == 2) && (std::string_view(argv[1]) == "--help"))
		{
			return 0;
		}

		return 1;
	}

	std::printf("Threading monkey test: duration=%us workers=%u max_threads=%u seed=%u\n",
	            state.options.duration, state.options.workers, state.options.max_threads, state.options.seed);
	std::fflush(stdout);

	// Include reference generation and worker shutdown in the watchdog's runtime limit
	state.workers = std::make_unique<WorkerStatus[]>(state.options.workers);
	std::mutex watchdog_lock;
	std::condition_variable watchdog_condition;
	bool is_finished = false;
	auto start_time = Clock::now();
	auto watchdog_deadline = start_time + std::chrono::seconds(state.options.duration + 30);
	std::thread watchdog([&] {
		std::unique_lock<std::mutex> guard(watchdog_lock);

		auto next_report = start_time + std::chrono::seconds(5);
		while (true)
		{
			auto wake_time = std::min(next_report, watchdog_deadline);
			bool is_completed = watchdog_condition.wait_until(guard, wake_time, [&] { return is_finished; });
			if (is_completed)
			{
				return;
			}

			auto now = Clock::now();
			if (now >= watchdog_deadline)
			{
				std::fprintf(stderr, "\nFAIL seed=%u: watchdog timeout\n",
				             state.options.seed);
				for (size_t i = 0; i < state.options.workers; i++)
				{
					std::fprintf(stderr, "Worker %zu iteration=%zu reference=%zu: %s\n", i,
					             state.workers[i].iteration.load(), state.workers[i].reference.load(),
					             state.workers[i].operation.load());
				}

				std::fflush(stderr);
				std::_Exit(1);
			}

			// Console output must not hold up the main thread's shutdown notification
			guard.unlock();
			print_statistics(state, std::chrono::duration<double>(now - start_time).count());
			guard.lock();
			next_report = Clock::now() + std::chrono::seconds(5);
		}
	});

	// Build immutable references before starting the randomized concurrent workloads
	make_images(state);
	if (make_references(state))
	{
		std::printf("Checking %zu image/settings pairs against serial references\n", state.references.size());
		std::fflush(stdout);
		auto deadline = Clock::now() + std::chrono::seconds(state.options.duration);
		std::vector<std::thread> workers;
		for (size_t i = 0; i < state.options.workers; i++)
		{
			workers.emplace_back(run_worker, std::ref(state), i, deadline);
		}

		for (auto& worker : workers)
		{
			worker.join();
		}
	}

	// Drain workers and release shared contexts before stopping the watchdog
	state.parents.clear();
	{
		std::lock_guard<std::mutex> guard(watchdog_lock);

		is_finished = true;
	}
	watchdog_condition.notify_one();
	watchdog.join();
	print_statistics(state, std::chrono::duration<double>(Clock::now() - start_time).count());
	print_references(state);

	if (state.is_failed)
	{
		return 1;
	}

	std::puts("\nPASS");
	return 0;
}
