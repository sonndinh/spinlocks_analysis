#include <iostream>
#include <cmath>
#include <random>
#include <map>
#include <algorithm>
#include <cassert>

#define _CXX_11_
#ifdef _CXX_11_
#include <chrono>
#else
#include <sys/time.h>
#endif
#include "TaskSetGenerator.h"

using namespace std;

//#define _DEBUG_DETERMINISTIC_
#ifdef _DEBUG_DETERMINISTIC_
unsigned int seed = 5;
unsigned int get_seed() {
	return seed++;
}

#else
unsigned int get_seed() {
#ifdef _CXX_11_
	unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
#else
	timeval time;
	gettimeofday(&time, NULL);
	unsigned int seed = time.tv_sec*1000000 + time.tv_usec;
#endif // _CXX_11_
	return seed;
}
#endif // _DEBUG_DETERMINISTIC_

/* Resource sharing factor */
double RSF[] = {0.1, 0.25, 0.4, 0.75};
//double RSF[] = {0.2, 0.25, 0.3, 0.4};


/**
 * Generate a task with period, critical path length, utilization
 * Use log-uniform distribution for period
 */
Task* create_task(unsigned int m) {
    mt19937 num_gen(get_seed());
	//	uniform_real_distribution<double> period_dist(MIN_TI, MAX_TI);
	uniform_real_distribution<double> log_period_dist(log10(MIN_TI), log10(MAX_TI));	

	Task *task = new Task();
	/* Generate period */
	//	double period = period_dist(num_gen);
	double log_period = log_period_dist(num_gen);
	double period = pow(10.0, log_period);
	if (ceil(period*1000) > MAX_TI*1000)
		task->T = floor(period*1000);
	else 
		task->T = ceil(period*1000);

	/* Generate critical path length */
	double critical_path_ratio;
	//	if (task->T >= (MAX_TI+MIN_TI)/2) {
	//		uniform_real_distribution<double> critical_path_dist(0.3, 0.425);
	//		critical_path_ratio = critical_path_dist(num_gen);
	//	} else {
	uniform_real_distribution<double> critical_path_dist(0.125, 0.25);
	critical_path_ratio = critical_path_dist(num_gen);
	//	}
	task->L = ceil(critical_path_ratio*task->T);

	/* Generate utilization */
	uniform_real_distribution<double> util_dist(1.25, sqrt(sqrt(m)));
	task->U = util_dist(num_gen);
	task->C = task->U * task->T;

	/* Calculate initial number of processors */
	task->procNum = ceil((task->C - task->L)/(task->T - task->L));
	
	//	cout << "Task parameters: <T,C,L,U> == <" << task->T << ", " << task->C 
	//		 << ", " << task->L << ", " << task->U << ">; #processors: " << task->procNum << endl;
	return task;
}


#if 0
Task* create_task(unsigned int m) {
	unsigned int seed = get_seed();
	/* C++0x does not have knuth_b generator, but C++11 does */
	knuth_b generator(seed);
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
	else
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
	
	//	cout << "Task parameters: <T,C,L,U> == <" << task->T << ", " << task->C 
	//		 << ", " << task->L << ", " << task->U << ">" << endl;
	//	cout << "Number of attempts for T: " << count1 << "; L: " << count2
	//		 << "; U:" << count3 << endl;
	return task;
}
#endif


/**
 * Generate taskset with bounded total utilization
 * Critical section lengths are generated with normal distributions
 */
#if 0
TaskSet* create_taskset(unsigned int m, unsigned int resourceNum, unsigned int Nmax, 
						CriticalDuration cslen_type) {
	TaskSet* tset = new TaskSet();
	double total_util = 0;

	/* Add new tasks until reach total utilization bound */
	int count_task = 0;
	const double MAX_UTIL = (double)m/2;
	while (total_util < MAX_UTIL) {
		Task* task = create_task(m);
		if ((total_util + task->U) > MAX_UTIL) {
			free(task);
			break;
		}
		
		count_task++;
		task->taskID = count_task;
		tset->tasks.insert(std::pair<TaskID, Task*>(task->taskID, task));
		total_util += task->U;
	}

	cout << "Number of tasks: " << tset->tasks.size() << endl;
	cout << "Total util: " << total_util << ". Max util: " << MAX_UTIL << endl;

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
		if (dirty_task_num <= 1) dirty_task_num = 2;
		//		cout << "Number of dirty tasks for resource " << i << ": " << dirty_task_num << endl;

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
			//			cout << "Task ID: " << dirty_task_id << " == <ResourceID: " << i << ", CSlen: " << res->CSLength 
			//				 << ", ReqNum: " << res->requestNum << ">" << endl;
		}
	}

	return tset;
}
#endif


/* Order tasks by decreasing periods (implicit deadline) */
bool compare(Task* first, Task* second) {
	return (first->T > second->T);
}

/**
 * Generate taskset with bounded total processors allocated
 * Critical section lengths are generated with uniform distributions
 */
TaskSet* create_taskset(unsigned int m, unsigned int resourceNum, unsigned int Nmax, 
						CriticalDuration cslen_type) {
	TaskSet* tset = new TaskSet();

	/* Add new tasks until the number of allocated processors reaches the bound */
	unsigned int count_task = 0;
	unsigned int proc_allocated = 0;
	//	const unsigned int MAX_PROC = ceil(m*0.95);
	const unsigned int MAX_PROC = m;
	while (proc_allocated < MAX_PROC) {
		Task * task = create_task(m);
		if (proc_allocated + task->procNum > MAX_PROC) {
			free(task);
			break;
		}

		count_task++;
		task->taskID = count_task;
		tset->tasks.insert(std::pair<TaskID, Task*> (task->taskID, task));
		proc_allocated += task->procNum;
	}

	//	cout << "Number of tasks: " << tset->tasks.size() << endl;

	unsigned int rsf_num = sizeof(RSF)/sizeof(RSF[0]);
	unsigned int task_num = tset->tasks.size();

	/* Set type of critical section length */
	double min_cslen, max_cslen;
	if (cslen_type == SHORT) {
		min_cslen = MIN_SHORT_CSLEN;
		max_cslen = MAX_SHORT_CSLEN;
	} else if (cslen_type == MODERATE) { 
		min_cslen = MIN_MODERATE_CSLEN;
		max_cslen = MAX_MODERATE_CSLEN;
	} else {
		min_cslen = MIN_LONG_CSLEN;
		max_cslen = MAX_LONG_CSLEN;
	}

	/* Uniform distribution of critical section lengths */
	mt19937 rangen(get_seed());
	//	minstd_ran rangen(get_seed());
	uniform_int_distribution<int> cslen_dist(min_cslen, max_cslen);

	/* Each resource has a resource-sharing-factor chosen randomly from RSF array */
	for (unsigned int i=1; i<=resourceNum; i++) {
		unsigned int rsf_idx = rangen() % rsf_num;
		unsigned int dirty_task_num = (unsigned int)(RSF[rsf_idx]*task_num);
		if (dirty_task_num <= 1)
			dirty_task_num = 2;

		//		cout << "Rsf production: " << RSF[rsf_idx]*task_num << "; Dirty task num: " << dirty_task_num << endl;
		//		cout << "Number of dirty tasks for resource " << i << ": " << dirty_task_num << endl;

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

			/* Generate number of accesses and critical section length */
			Resource* res = new Resource();
			res->resourceID = i;
			res->requestNum = rangen() % Nmax + 1;
			res->CSLength = cslen_dist(rangen);

			/* Add to resource list of the task */
			tset->tasks[dirty_task_id]->myResources.insert(std::pair<ResourceID, Resource*>(i, res));

			//			cout << "Task ID: " << dirty_task_id << " == <ResourceID: " << i << ", CSlen: " << res->CSLength 
			//				 << ", ReqNum: " << res->requestNum << ">" << endl;
		}
	}


	/* Assign tasks' priorities (a larger number means a higher priority) */
	vector<Task*> tasks_array;
	for (map<TaskID, Task*>::iterator it=tset->tasks.begin(); it!=tset->tasks.end(); it++) {
		tasks_array.push_back(it->second);
	}

	/* Sort tasks by decreasing relative deadline */
	sort(tasks_array.begin(), tasks_array.end(), compare);
	unsigned int prio = 1;
	tasks_array[0]->priority = prio;
	for (int i=1; i<tasks_array.size(); i++) {
		if (tasks_array[i]->T == tasks_array[i-1]->T) {
			tasks_array[i]->priority = tasks_array[i-1]->priority;
		} else {
			tasks_array[i]->priority = tasks_array[i-1]->priority + 1;
		}
	}


	/* Assert priority assignment */
	/*
	cout << "Task ID: " << tasks_array[0]->taskID << "; Deadline: " << tasks_array[0]->T 
		 << "; Priority: " << tasks_array[0]->priority << endl;
	for (int i=1; i<tasks_array.size(); i++) {
		assert(tasks_array[i]->T <= tasks_array[i-1]->T);
		assert(tasks_array[i]->priority >= tasks_array[i-1]->priority);
		cout << "Task ID: " << tasks_array[i]->taskID << "; Deadline: " << tasks_array[i]->T 
			 << "; Priority: " << tasks_array[i]->priority << endl;
	}
	*/

	/* For each task, find the sets of higher, lower, and equal priority tasks */
	for (map<TaskID, Task*>::iterator it=tset->tasks.begin(); it!=tset->tasks.end(); it++) {
		TaskID task_id = it->first;
		vector<TaskID> &hp_tasks = tset->higher_prio_tasks[task_id];
		vector<TaskID> &lp_tasks = tset->lower_prio_tasks[task_id];
		vector<TaskID> &ep_tasks = tset->equal_prio_tasks[task_id];

		for (int i=0; i<tasks_array.size(); i++) {
			if (tasks_array[i]->taskID == task_id)
				continue;

			if (tasks_array[i]->priority > it->second->priority) {
				hp_tasks.push_back(tasks_array[i]->taskID);
			} else if(tasks_array[i]->priority < it->second->priority) {
				lp_tasks.push_back(tasks_array[i]->taskID);
			} else {
				ep_tasks.push_back(tasks_array[i]->taskID);
			}
		}
	}

	return tset;
}
