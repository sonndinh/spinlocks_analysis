#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <algorithm>

#define _CXX_11_
#ifdef _CXX_11_
#include <chrono>
#else
#include <sys/time.h>
#endif
#include "TaskSetGenerator.h"

using namespace std;

unsigned int seed = 5;
unsigned int get_seed() {
	return seed++;
}

#if 0
unsigned int get_seed() {
#ifdef _CXX_11_
	unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
#else
	timeval time;
	gettimeofday(&time, NULL);
	unsigned int seed = time.tv_sec*1000000 + time.tv_usec;
#endif
	return seed;
}
#endif

/* Resource sharing factor */
double RSF[] = {0.2, 0.3, 0.5, 0.75};

Task* create_task(unsigned int m) {
	unsigned int seed = get_seed();
	/* C++0x does not have knuth_b generator, but C++11 does */
	knuth_b generator(seed);
	//	minstd_rand generator(seed);
	normal_distribution<double> period_dist(MEAN_TI, STDEV_TI);

	int count1,count2,count3;
	count1 = count2 = count3 = 1;
	Task* task = new Task();

	/* Generate task period */
	double t_i = period_dist(generator);
	while(t_i < MIN_TI || t_i > MAX_TI) {
		// Get new value for period
		t_i = period_dist(generator);
		count1++;
	}
	if (ceil(t_i*1000) > MAX_TI*1000) 
		task->T = floor(t_i*1000);
	task->T = ceil(t_i*1000);

	/* Generate critical path length parameters:
	 * max, min, mean, stdev (ms)
	 */
	double max_li = t_i/3;
	double min_li = 1;
	double mean_li = (min_li + max_li)/2;
	double stdev_li = (min_li + max_li)/4;
	normal_distribution<double> criticalpath_dist(mean_li, stdev_li);
	double l_i = criticalpath_dist(generator);
	while (l_i < min_li || l_i > max_li) {
		l_i = criticalpath_dist(generator);
		count2++;
	}
	task->L = ceil(l_i*1000);

	/* Generate task utilization in range (1; sqrt(m)/3) 
	 */
	double min_ui = 1;
	/* Be careful: if sqrt(m) <= 3, then max_ui is <= 1 */
	double max_ui = sqrt(m)/3;
	double mean_ui = (min_ui + max_ui)/2;
	double stdev_ui = (min_ui + max_ui)/4;
	normal_distribution<double> utilization_dist(mean_ui, stdev_ui);
	double u_i = utilization_dist(generator);
	while (u_i <= min_ui || u_i >= max_ui) {
		u_i = utilization_dist(generator);
		count3++;
	}
	task->U = u_i;
	task->C = u_i*task->T;
	
	cout << "Task parameters: <T,C,L,U> == <" << task->T << ", " << task->C 
		 << ", " << task->L << ", " << task->U << ">" << endl;
	//	cout << "Number of attempts for T: " << count1 << "; L: " << count2
	//		 << "; U:" << count3 << endl;
	return task;
}


TaskSet* create_taskset(unsigned int m, unsigned int resourceNum, unsigned int Nmax, 
						CriticalDuration cslen_type) {
	TaskSet* tset = new TaskSet();
	double total_util = 0;

	/* Add new tasks until reach total utilization bound */
	int count_task = 1;
	while (total_util < (double)m/2) {
		Task* task = create_task(m);
		if ((total_util + task->U) > (double)m/2) {
			free(task);
			break;
		}
		
		task->taskID = count_task;
		count_task++;
		tset->tasks.insert(std::pair<TaskID, Task*>(task->taskID, task));
		total_util += task->U;
	}

	cout << "Number of tasks: " << tset->tasks.size() << endl;
	cout << "Total util: " << total_util << ". Max util: " << m/2 << endl;

	unsigned int rsf_num = sizeof(RSF)/sizeof(RSF[0]);
	unsigned int task_num = tset->tasks.size();

	/* Common random number generator for choosing rsf, dirty task */
	unsigned int seed = get_seed();
	minstd_rand rangen(seed);
	/* Set type of critical section length */
	double mean_cslen, stdev_cslen, min_cslen, max_cslen;
	if (cslen_type == Short) {
		mean_cslen = MEAN_SHORT_CSLEN;
		stdev_cslen = STDEV_SHORT_CSLEN;
		min_cslen = MIN_SHORT_CSLEN;
		max_cslen = MAX_SHORT_CSLEN;
	} else { 
		mean_cslen = MEAN_LONG_CSLEN;
		stdev_cslen = STDEV_LONG_CSLEN;
		min_cslen = MIN_LONG_CSLEN;
		max_cslen = MAX_LONG_CSLEN;
	}
	normal_distribution<double> csLen_dist(mean_cslen, stdev_cslen);

	/* Each resource has a rsf chosen randomly from RSF array */
	for (unsigned int i=1; i<=resourceNum; i++) {
		unsigned int rsf_idx = rangen() % rsf_num;
		unsigned int dirty_task_num = (unsigned int)(RSF[rsf_idx]*task_num);
		if (dirty_task_num == 1) dirty_task_num = 2;
		cout << "Number of dirty tasks for resource " << i << ": " << dirty_task_num << endl;

		vector<TaskID> dirty_tasks;
		unsigned int count = 0;
		while (count < dirty_task_num) {
			/* Pick a random task from task set */
			unsigned int dirty_task_id = rangen() % task_num + 1;
			if (find(dirty_tasks.begin(), dirty_tasks.end(), dirty_task_id) != dirty_tasks.end()) {
				continue;
			} else {
				count++;
				dirty_tasks.push_back(dirty_task_id);
			}
			
			//cout << "Picked task: " << dirty_task_id << endl;
			Resource* res = new Resource();
			res->resourceID = i;
			res->requestNum = rangen() % Nmax + 1;
			res->CSLength = ceil(csLen_dist(rangen));
			while (res->CSLength < min_cslen || res->CSLength > max_cslen) {
				res->CSLength = ceil(csLen_dist(rangen));
			}
			
			tset->tasks[dirty_task_id]->myResources.insert(std::pair<ResourceID, Resource*>(i, res));
			cout << "Task ID: " << dirty_task_id << " == <ResourceID: " << i << ", CSlen: " << res->CSLength 
				 << ", ReqNum: " << res->requestNum << ">" << endl;
		}
	}

	return tset;
}
