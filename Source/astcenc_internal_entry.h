// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2011-2026 Arm Limited
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
 * @brief Functions and data declarations for the outer context.
 *
 * The outer context includes thread-pool management, which is slower to
 * compile due to increased use of C++ stdlib. The inner context used in the
 * majority of the codec library does not include this.
 */

#if !defined(ASTCENC_INTERNAL_ENTRY_INCLUDED)
#define ASTCENC_INTERNAL_ENTRY_INCLUDED

#include <condition_variable>
#include <functional>
#include <mutex>

#include "astcenc_internal.h"

/* ============================================================================
  Parallel execution control
============================================================================ */

/**
 * @brief A simple counter-based manager for parallel task execution.
 *
 * The task processing execution consists of:
 *
 *     * A single-threaded init stage.
 *     * A multi-threaded processing stage.
 *     * A condition variable so threads can wait for processing completion.
 *
 * The init stage will be executed by the first thread to arrive in the critical section, there is
 * no main thread in the thread pool.
 *
 * The processing stage uses dynamic dispatch to assign task tickets to threads on an on-demand
 * basis. Threads may each therefore executed different numbers of tasks, depending on their
 * processing complexity. The task queue and the task tickets are just counters; the caller must map
 * these integers to an actual processing partition in a specific problem domain.
 *
 * The exit wait condition is needed to ensure processing has finished before a worker thread can
 * progress to the next stage of the pipeline. Specifically a worker may exit the processing stage
 * because there are no new tasks to assign to it while other worker threads are still processing.
 * Calling @c wait() will ensure that all other worker have finished before the thread can proceed.
 *
 * The basic usage model:
 *
 *     // --------- From single-threaded code ---------
 *
 *     // Reset the tracker state
 *     manager->reset()
 *
 *     // --------- From multi-threaded code ---------
 *
 *     // Run the stage init; only first thread actually runs the lambda
 *     manager->init(<lambda>)
 *
 *     do
 *     {
 *         // Request a task assignment
 *         uint task_count;
 *         uint base_index = manager->get_tasks(<granule>, task_count);
 *
 *         // Process any tasks we were given (task_count <= granule size)
 *         if (task_count)
 *         {
 *             // Run the user task processing code for N tasks here
 *             ...
 *
 *             // Flag these tasks as complete
 *             manager->complete_tasks(task_count);
 *         }
 *     } while (task_count);
 *
 *     // Wait for all threads to complete tasks before progressing
 *     manager->wait()
 *
  *     // Run the stage term; only first thread actually runs the lambda
 *     manager->term(<lambda>)
 */
class ParallelManager
{
private:
	/** @brief Lock used for critical section and condition synchronization. */
	std::mutex m_lock;

	/** @brief True if the current operation is cancelled. */
	bool m_is_cancelled;

	/** @brief True if the stage init() step has been executed. */
	bool m_init_done;

	/** @brief True if the stage term() step has been executed. */
	bool m_term_done;

	/** @brief Condition variable for tracking stage processing completion. */
	std::condition_variable m_complete;

	/** @brief Number of tasks issued to workers, but not necessarily finished. */
	size_t m_issued_task_count;

	/** @brief Number of tasks finished. */
	size_t m_done_task_count;

	/** @brief Total number of tasks in the stage. */
	size_t m_total_task_count;

	/** @brief Progress callback (optional). */
	astcenc_progress_callback m_callback;

	/** @brief Lock used for callback synchronization. */
	std::mutex m_callback_lock;

	/** @brief Minimum progress before making a callback. */
	float m_callback_min_diff;

	/** @brief Last progress callback value, protected by @c m_callback_lock. */
	float m_callback_last_value;

public:
	/** @brief Create a new ParallelManager. */
	ParallelManager()
	{
		reset();
	}

	/**
	 * @brief Reset the tracker for a new processing batch.
	 *
	 * This must be called from single-threaded code before starting the multi-threaded processing
	 * operations.
	 */
	void reset()
	{
		// Single-threaded contexts reset on entry while another thread may request cancellation
		std::lock_guard<std::mutex> lck(m_lock);

		m_init_done = false;
		m_term_done = false;
		m_is_cancelled = false;
		m_issued_task_count = 0;
		m_done_task_count = 0;
		m_total_task_count = 0;
		m_callback = nullptr;
		m_callback_last_value = 0.0f;
		m_callback_min_diff = 1.0f;
	}

	/**
	 * @brief Stop new tasks being assigned while preserving outstanding work.
	 *
	 * Note, all in-flight tasks in a worker will still complete normally.
	 */
	void cancel()
	{
		std::unique_lock<std::mutex> lck(m_lock);

		m_is_cancelled = true;
		lck.unlock();
		m_complete.notify_all();
	}

	/**
	 * @brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to hit this will process the
	 * initialization. Other threads will block and wait for it to complete.
	 *
	 * @param init_func   Callable which executes the stage initialization. It must return the
	 *                    total number of tasks in the stage.
	 */
	void init(std::function<size_t(void)> init_func)
	{
		std::lock_guard<std::mutex> lck(m_lock);

		if (!m_init_done)
		{
			m_total_task_count = init_func();
			m_init_done = true;
		}
	}

	/**
	 * @brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to hit this will process the
	 * initialization. Other threads will block and wait for it to complete.
	 *
	 * @param task_count   Total number of tasks needing processing.
	 * @param callback     Function pointer for progress status callbacks.
	 */
	void init(size_t task_count, astcenc_progress_callback callback)
	{
		std::lock_guard<std::mutex> lck(m_lock);

		if (!m_init_done)
		{
			m_callback = callback;
			m_total_task_count = task_count;
			m_init_done = true;

			// Report every 1% or 4096 blocks, whichever is larger, to avoid callback overhead
			float min_diff = (4096.0f / static_cast<float>(task_count)) * 100.0f;
			m_callback_min_diff = astc::max(min_diff, 1.0f);
		}
	}

	/**
	 * @brief Request a task assignment.
	 *
	 * Assign up to @c granule tasks to the caller for processing.
	 *
	 * @param      granule   Maximum number of tasks that can be assigned.
	 * @param[out] count     Actual number of tasks assigned, or zero if no tasks were assigned.
	 *
	 * @return Task index of the first assigned task; assigned tasks increment from this.
	 */
	size_t get_task_assignment(size_t granule, size_t& count)
	{
		std::lock_guard<std::mutex> lck(m_lock);

		if (m_is_cancelled || (m_issued_task_count >= m_total_task_count))
		{
			count = 0;
			return 0;
		}

		size_t base = m_issued_task_count;
		count = astc::min(m_total_task_count - base, granule);
		m_issued_task_count += count;
		return base;
	}

	/**
	 * @brief Complete a task assignment.
	 *
	 * Mark @c count tasks as complete. This will notify all threads blocked on @c wait() if this
	 * completes the processing of the stage.
	 *
	 * @param count   The number of completed tasks.
	 */
	void complete_task_assignment(size_t count)
	{
		// Update m_done_task_count under the mutex to synchronize with threads entering wait()
		size_t local_count;
		{
			std::unique_lock<std::mutex> lck(m_lock);
			m_done_task_count += count;
			local_count = m_done_task_count;

			// Notify if nothing left to do
			size_t task_count = m_total_task_count;
			if (m_is_cancelled)
			{
				task_count = m_issued_task_count;
			}

			if (m_done_task_count == task_count)
			{
				lck.unlock();
				m_complete.notify_all();
			}
		}

		// Report progress when a callback is registered
		if (m_callback)
		{
			// Keep callbacks outside m_lock so a callback can cancel compression
			std::lock_guard<std::mutex> cblck(m_callback_lock);

			float num = static_cast<float>(local_count);
			float den = static_cast<float>(m_total_task_count);
			float this_value =  (num / den) * 100.0f;
			bool report_test = (local_count == m_total_task_count) ||
			                   ((this_value - m_callback_last_value) > m_callback_min_diff);

			if (report_test)
			{
				m_callback(this_value);
				m_callback_last_value = this_value;
			}
		}
	}

	/**
	 * @brief Wait for issued tasks to finish, including after cancellation.
	 */
	void wait()
	{
		std::unique_lock<std::mutex> lck(m_lock);

		m_complete.wait(lck, [this] {
			size_t task_count = m_total_task_count;
			if (m_is_cancelled)
			{
				task_count = m_issued_task_count;
			}

			return m_done_task_count == task_count;
		});
	}

	/**
	 * @brief Trigger the pipeline stage term step.
	 *
	 * This can be called from multi-threaded code. The first thread to hit this will process the
	 * work pool termination. Caller must have called @c wait() prior to calling this function to
	 * ensure that processing is complete.
	 *
	 * @param term_func   Callable which executes the stage termination.
	 */
	void term(std::function<void(void)> term_func)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_term_done)
		{
			term_func();
			m_term_done = true;
		}
	}
};

/**
 * @brief The astcenc compression context.
 */
struct astcenc_context
{
	/** @brief The context internal state. */
	astcenc_contexti context;

#if !defined(ASTCENC_DECOMPRESS_ONLY)
	/** @brief The parallel manager for averages computation. */
	ParallelManager manage_avg;

	/** @brief The parallel manager for compression. */
	ParallelManager manage_compress;
#endif

	/** @brief The parallel manager for decompression. */
	ParallelManager manage_decompress;
};

#endif
