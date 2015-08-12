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

#include "taskset_generator.h"
#include "blocking_analysis.h"

using namespace std;

// #define _DEBUG_DETERMINISTIC_
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

// Resource sharing factor
double g_rsf[] = {0.1, 0.25, 0.4, 0.75};
//double RSF[] = {0.2, 0.25, 0.3, 0.4};


// Generate a task with period, critical path length, utilization
// Use log-uniform distribution for period
Task* create_task(unsigned int m) {
    mt19937 num_gen(get_seed());
	uniform_real_distribution<double> log_period_dist(log10(MIN_TI), log10(MAX_TI));
	Task *task = new Task();

	// Generate period
	double log_period = log_period_dist(num_gen);
	double period = pow(10.0, log_period);
	if (ceil(period*1000) > MAX_TI*1000)
		task->period_ = floor(period*1000);
	else 
		task->period_ = ceil(period*1000);

	// Generate critical path length
	double critical_path_ratio;
	//	if (task->T >= (MAX_TI+MIN_TI)/2) {
	//		uniform_real_distribution<double> critical_path_dist(0.3, 0.425);
	//		critical_path_ratio = critical_path_dist(num_gen);
	//	} else {
	uniform_real_distribution<double> critical_path_dist(0.125, 0.25);
	critical_path_ratio = critical_path_dist(num_gen);
	//	}
	task->span_ = ceil(critical_path_ratio*task->period_);

	// Generate utilization
	uniform_real_distribution<double> util_dist(1.25, sqrt(sqrt(m)));
	task->utilization_ = util_dist(num_gen);
	task->wcet_ = task->utilization_ * task->period_;

	// Calculate initial number of processors
	task->proc_num_ = ceil((task->wcet_ - task->span_)/(task->period_ - task->span_));
	
	return task;
}


// Order tasks by decreasing periods (implicit deadline)
bool compare_period(const Task* first, const Task* second) {
	return (first->period_ > second->period_);
}

// Generate taskset with bounded total processors allocated
// Critical section lengths are generated with uniform distributions
TaskSet* create_taskset(unsigned int m, unsigned int resource_num, unsigned int n_max, 
						CriticalDuration cslen_type) {
	TaskSet* tset = new TaskSet();

	// Add new tasks until the number of allocated processors reaches the bound
	unsigned int count_task = 0;
	unsigned int proc_allocated = 0;
	//	const unsigned int MAX_PROC = ceil(m*0.95);
	const unsigned int MAX_PROC = m;
	while (proc_allocated < MAX_PROC) {
		Task * task = create_task(m);
		if (proc_allocated + task->proc_num_ > MAX_PROC) {
			delete task;
			break;
		}

		count_task++;
		task->task_id_ = count_task;
		tset->tasks_.insert(std::pair<TaskID, Task*> (task->task_id_, task));
		proc_allocated += task->proc_num_;
	}

	//	cout << "Number of tasks: " << tset->tasks.size() << endl;

	unsigned int rsf_num = sizeof(g_rsf)/sizeof(g_rsf[0]);
	unsigned int task_num = tset->tasks_.size();

	// Set type of critical section length
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

	// Uniform distribution of critical section lengths
	mt19937 rangen(get_seed());
	uniform_int_distribution<int> cslen_dist(min_cslen, max_cslen);

	// Each resource has a resource-sharing-factor chosen randomly from RSF array
	for (unsigned int i=1; i<=resource_num; i++) {
		unsigned int rsf_idx = rangen() % rsf_num;
		unsigned int dirty_task_num = (unsigned int)(g_rsf[rsf_idx]*task_num);
		if (dirty_task_num <= 1)
			dirty_task_num = 2;

		//		cout << "Rsf production: " << RSF[rsf_idx]*task_num << "; Dirty task num: " << dirty_task_num << endl;
		//		cout << "Number of dirty tasks for resource " << i << ": " << dirty_task_num << endl;

		vector<TaskID> dirty_tasks;
		unsigned int count = 0;
		while (count < dirty_task_num) {
			// Pick a random task from task set
			unsigned int dirty_task_id = rangen() % task_num + 1;
			if (find(dirty_tasks.begin(), dirty_tasks.end(), dirty_task_id) != dirty_tasks.end()) {
				continue;
			} else {
				count++;
				dirty_tasks.push_back(dirty_task_id);
			}
			
			//cout << "Picked task: " << dirty_task_id << endl;

			// Generate number of accesses and critical section length
			Resource* res = new Resource();
			res->resource_id = i;
			res->request_num = rangen() % n_max + 1;
			res->critical_section_len = cslen_dist(rangen);

			// Add to resource list of the task
			tset->tasks_[dirty_task_id]->my_resources_.insert(std::pair<ResourceID, Resource*>(i, res));
		}
	}


	// Assign tasks' priorities (a larger number means a higher priority)
	vector<Task*> tasks_array;
	for (map<TaskID, Task*>::iterator it=tset->tasks_.begin(); it!=tset->tasks_.end(); it++) {
		tasks_array.push_back(it->second);
	}

	// Sort tasks by decreasing relative deadline
	sort(tasks_array.begin(), tasks_array.end(), compare_period);
	unsigned int priority = 1;
	tasks_array[0]->priority_ = priority;
	for (int i=1; i<tasks_array.size(); i++) {
		if (tasks_array[i]->period_ == tasks_array[i-1]->period_) {
			tasks_array[i]->priority_ = tasks_array[i-1]->priority_;
		} else {
			tasks_array[i]->priority_ = tasks_array[i-1]->priority_ + 1;
		}
	}


	// Assert priority assignment
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

	// For each task, find the sets of higher, lower, and equal priority tasks
	for (map<TaskID, Task*>::iterator it=tset->tasks_.begin(); it!=tset->tasks_.end(); it++) {
		TaskID task_id = it->first;
		vector<TaskID> &hp_tasks = tset->higher_prio_tasks_[task_id];
		vector<TaskID> &lp_tasks = tset->lower_prio_tasks_[task_id];
		vector<TaskID> &ep_tasks = tset->equal_prio_tasks_[task_id];

		for (int i=0; i<tasks_array.size(); i++) {
			if (tasks_array[i]->task_id_ == task_id)
				continue;

			if (tasks_array[i]->priority_ > it->second->priority_) {
				hp_tasks.push_back(tasks_array[i]->task_id_);
			} else if(tasks_array[i]->priority_ < it->second->priority_) {
				lp_tasks.push_back(tasks_array[i]->task_id_);
			} else {
				ep_tasks.push_back(tasks_array[i]->task_id_);
			}
		}
	}

	return tset;
}

// Initialize blocking time, #processors, response time 
// Also, update "interferences" for each task in task set
// Return: true, if taskset are OK after initiating
//         false, otherwise
// (OK means total processors allocated to tasks < m AND 
//  there is no task with Response time > Deadline)
bool init_iteration(TaskSet* taskset, unsigned int m) {
	map<TaskID, Task*> &tset = taskset->tasks_;
	map<TaskID, Task*>::iterator it = tset.begin();
	for (; it != tset.end(); it++) {
		Task* task = it->second;
		task->blocking_on_single_ = 0;
		task->proc_num_ = BlockingAnalysis::get_instance().alloc_proc(task->wcet_, task->span_, task->period_, 0, 0);
		task->response_time_ = BlockingAnalysis::get_instance().response_time(task->wcet_, task->span_, task->proc_num_, 0, 0);
		task->converged_ = false;
		
		// Update interferences for this task
		map<ResourceID, vector<TaskID>*> &interferences = task->interferences_;
		map<ResourceID, Resource*> &my_resources = task->my_resources_;
		map<TaskID, Task*>::iterator task_iter;
		map<ResourceID, Resource*>::iterator res_iter;

		for (res_iter = my_resources.begin(); res_iter!=my_resources.end(); res_iter++) {
			ResourceID rid = res_iter->first;
			interferences.insert(std::pair<ResourceID, vector<TaskID>*> (rid, new vector<TaskID>()));
			for (task_iter = tset.begin(); task_iter!=tset.end(); task_iter++) {
				// Abort if this is me
				if (task_iter->second->task_id_ == task->task_id_)
					continue;

				// Add tasks access to the same resource
				if (task_iter->second->my_resources_.find(rid) != 
					task_iter->second->my_resources_.end()) {
					interferences.find(rid)->second->push_back(task_iter->first);
				}
			}
		}
	}

	/*
	//#define _INIT_DEBUG_
	//#define _INIT_VERBOSE_
#ifdef _INIT_DEBUG_
	// Debug: dump tasks information
	unsigned int total_proc = 0;
	bool isOk = true;
	int number_tasks_fail = 0;
	for (it=tset.begin(); it!=tset.end(); it++) {
#ifdef _INIT_VERBOSE_
		cout << "Task ID: " << it->second->taskID << endl;
		cout << "Parameters (T,L,U,C): (" << it->second->T << "," <<
			it->second->L << "," << it->second->U << "," << 
			it->second->C << ")" << endl;
		cout << "#processors: " << it->second->procNum << endl;
		cout << "Response time: " << it->second->R << endl;
		cout << "Response time > Deadline: " <<
			((it->second->R > it->second->T)? "true" : "false") << endl;
#endif // _INIT_VERBOSE_
		
		total_proc += it->second->procNum;
		if (it->second->R > it->second->T) {
			number_tasks_fail++;
			isOk = false;
		}
		map<ResourceID, Resource*> &res = it->second->myResources;
		map<ResourceID, Resource*>::iterator resIt = res.begin();
#ifdef _INIT_VERBOSE_
		cout << "Resource ID list: (";
		for (; resIt!=res.end(); resIt++) {
			cout << resIt->first << " ";
		}
		cout << ")" << endl;

		map<ResourceID, vector<TaskID>*> &interfe = it->second->interferences;
		map<ResourceID, vector<TaskID>*>::iterator intIt = interfe.begin();
		cout << "Interference List ============= " << endl;
		for (; intIt!=interfe.end(); intIt++) {
			cout << "   For resource ID " << intIt->first << ": ";
			vector<TaskID>* vec = intIt->second;
			for (int i=0; i<vec->size(); i++) {
				cout << vec->at(i) << " ";
			}
			cout << endl;
		}
		cout << endl;
#endif // _INIT_VERBOSE_
	}
	cout << "Total #processors (initialization): " << total_proc << endl;
	if (total_proc > m || isOk == false) {
		cout << "#tasks fail: " << number_tasks_fail << "; Taskset fail !!!" << endl;
		return false;
	}
	return true;
#endif // _INIT_DEBUG_
	*/

	// In case of no debugging, it just returns true 
	// But it does not necessarily mean the task set is OK
	return true;
}
