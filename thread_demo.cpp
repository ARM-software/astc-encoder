// SPDX-License-Identifier: Apache-2.0
// ----------------------------------------------------------------------------
// Copyright 2020 Arm Limited
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
 * @brief Prototype worker thread manager.
 */

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <iostream>
#include <mutex>
#include <thread>

/**
 * @brief A simple counter-based manager for parallel task execution.
 *
 * The task processing execution consists of:
 *
 *     * A single-threaded init stage.
 *     * A multi-threaded processing stage.
 *     * A condition variable so threads can wait for processing completion.
 *
 * The init stage will be executed by the first thread to arrive in the
 * critical section, there is no master thread in the thread pool.
 *
 * The processing stage uses dynamic dispatch to assign task tickets to threads
 * on an on-demand basis. Threads may each therefore executed different numbers
 * of tasks, depending on their processing complexity. The task queue and the
 * task tickets are just counters; the caller must map these integers to an
 * actual processing partition in a specific problem domain.
 *
 * The exit wait condition is needed to ensure processing has finished before
 * a worker thread can progress to the next stage of the pipeline. Specifically
 * a worker may exit the processing stage because there are no new tasks to
 * assign to it while other worker threads are still processing. Calling wait()
 * will ensure that all other worker have finished before the thread can
 * proceed.
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
 *     manager->init(<total_task_count>, <lambda>)
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
 *             // Run the user task processing code here
 *             ...
 *
 *             // Flag these tasks as complete
 *             manager->complete_tasks(task_count);
 *         }
 *     } while (task_count);
 *
 *     // Wait for all threads to complete tasks before progressing
 *     manager->wait()
 */
class ParallelManager
{
private:
	/** \brief Lock used for critical section and condition synchronization. */
	std::mutex m_lock;

	/** \brief True if the stage init() step has been executed. */
	bool m_init_done;

	/** \brief Contition variable for tracking stage processing completion. */
	std::condition_variable m_complete;

	/** \brief Number of tasks started, but not necessarily finished. */
	unsigned int m_start_count;

	/** \brief Number of tasks finished. */
	unsigned int m_done_count;

	/** \brief Number of tasks that need to be processed. */
	unsigned int m_task_count;

public:
	/** \brief Create a new ParallelManager. */
	ParallelManager()
	{
		reset();
	}

	/**
	 * \brief Reset the tracker for a new processing batch.
	 *
	 * This must be called from single-threaded code before starting the
	 * multi-threaded procesing operations.
	 */
	void reset()
	{
		m_init_done = false;
		m_start_count = 0;
		m_done_count = 0;
		m_task_count = 0;
	}

	/**
	 * \brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to
	 * hit this will process the initialization. Other threads will block
	 * and wait for it to complete.
	 *
	 * \param task_count   Total number of tasks needing processing.
	 * \param init_func    Callable which executes the stage initialization.
	 */
	void init(unsigned int task_count, std::function<void(void)> init_func)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_init_done)
		{
			init_func();
			m_task_count = task_count;
			m_init_done = true;
		}
	}

	/**
	 * \brief Trigger the pipeline stage init step.
	 *
	 * This can be called from multi-threaded code. The first thread to
	 * hit this will process the initialization. Other threads will block
	 * and wait for it to complete.
	 *
	 * \param task_count   Total number of tasks needing processing.
	 */
	void init(unsigned int task_count)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		if (!m_init_done)
		{
			m_task_count = task_count;
			m_init_done = true;
		}
	}

	/**
	 * \brief Request a task assignment.
	 *
	 * Assign up to \c granule tasks to the caller for processing.
	 *
	 * \param      granule   Maximum number of tasks that can be assigned.
	 * \param[out] count     Actual number of tasks assigned, or zero if
	 *                       no tasks were assigned.
	 *
	 * \return Task index of the first assigned task; assigned tasks
	 *         increment from this.
	 */
	unsigned int get_task_assignment(unsigned int granule, unsigned int& count)
	{
		std::lock_guard<std::mutex> lck(m_lock);
		unsigned int base = m_start_count;
		count = std::min(granule, m_task_count - m_start_count);
		m_start_count += count;
		return base;
	}

	/**
	 * \brief Complete a task assignment.
	 *
	 * Mark \c count tasks as complete. This will notify all threads blocked
	 * on \c wait() if this completes the processing of the stage.
	 *
	 * \param count   The number of completed tasks.
	 */
	void complete_task_assignment(unsigned int count)
	{
		std::unique_lock<std::mutex> lck(m_lock);
		this->m_done_count += count;
		if (m_done_count == m_task_count)
		{
			lck.unlock();
			m_complete.notify_all();
		}
	}

	/**
	 * \brief Wait for stage processing to complete.
	 */
	void wait()
	{
		std::unique_lock<std::mutex> lck(m_lock);
		m_complete.wait(lck, [this]{ return m_done_count == m_task_count; });
	}
};

void run(int tid, unsigned int total_count, ParallelManager* manager)
{
	const char* foo = "Demo init";
	auto initfunc = [foo]() { std::cout << foo << "\n"; };

	unsigned int base_index;
	unsigned int task_count;
	manager->init(total_count, initfunc);

	do
	{
		base_index = manager->get_task_assignment(2, task_count);
		if (task_count)
		{
			// TODO - Processing would go here
			manager->complete_task_assignment(task_count);
		}
	} while (task_count);

	manager->wait();
}

int main(void)
{
	ParallelManager manager;
	unsigned int total_count = 21;

	std::thread th0(run, 0, total_count, &manager);
	std::thread th1(run, 1, total_count, &manager);
	std::thread th2(run, 2, total_count, &manager);

	th0.join();
	th1.join();
	th2.join();

	return 0;
}
