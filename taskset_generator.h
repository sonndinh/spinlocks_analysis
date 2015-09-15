#ifndef __TASKSET_GENERATOR_H__
#define __TASKSET_GENERATOR_H__

#include <map>
#include <vector>
#include "defines.h"

using namespace std;

// Structure for storing resource data under a particular task
typedef struct Resource {
	ResourceID resource_id;
	double critical_section_len; // critical section length, in microsecond
	unsigned int request_num;  // number of requests from an enclosed task
} Resource;

// Implicit deadline parallel task
class Task { 
 public:
	// Copy assignment operator which is used to clone a task set
	Task &operator= (const Task& other) {
		if (this != &other) {
			task_id_ = other.task_id_;
			period_ = other.period_;
			span_ = other.span_;
			utilization_ = other.utilization_;
			wcet_ = other.wcet_;
			
			blocking_on_single_ = other.blocking_on_single_;
			proc_num_ = other.proc_num_;
			response_time_ = other.response_time_;
			priority_ = other.priority_;
			
			map<ResourceID, Resource*>::const_iterator res_it;
			for (res_it = other.my_resources_.begin(); res_it != other.my_resources_.end(); res_it++) {
				Resource *my_res = new Resource();
				*my_res = *(res_it->second);
				my_resources_.insert( std::pair<ResourceID, Resource*>(res_it->first, my_res) );
			}

			map<ResourceID, vector<TaskID>*>::const_iterator inter_it;
			for (inter_it = other.interferences_.begin(); inter_it!=other.interferences_.end(); inter_it++) {
				vector<TaskID>* my_vec = new vector<TaskID>;
				*my_vec = *(inter_it->second);
				interferences_.insert( std::pair<ResourceID, vector<TaskID>*>(inter_it->first, my_vec) );
			}
		}
		
		return *this;
	};

	// Destructor
	~Task() {
		for (map<ResourceID, Resource*>::iterator it=my_resources_.begin(); it!=my_resources_.end(); it++) {
			delete it->second;
		}

		for (map<ResourceID, vector<TaskID>*>::iterator it=interferences_.begin(); it!=interferences_.end(); it++) {
			delete it->second;
		}
	}

 public:
	TaskID task_id_;
	double period_; // period, in microsecond
	double span_; // critical path length, in microsecond
	double utilization_; // task's utiliation
	double wcet_; // Worst-case-execution-time
	
	map<ResourceID, Resource*> my_resources_; // map resouce id to its data

	// Stored data for fixed-point iteration
	double blocking_on_single_; // worst-case blocking time on 1 processor, in microsecond
	unsigned int proc_num_; // number of processors allocated
	double response_time_; // response time, in microsecond

	unsigned int new_proc_num_; // new number of processors allocated in the next iteration
	bool converged_; // true if the task parameters have converged
	
	// Other tasks access to the same resources
	map<ResourceID, vector<TaskID>*> interferences_;

	// Larger number means higher priority
	unsigned int priority_;
};


// Task set comprised of parallel tasks
class TaskSet {
 public:
	// Copy assignment operator to clone a task set
	TaskSet &operator= (const TaskSet &other) {
		if (this != &other) {
			// Copy list of tasks
			map<TaskID, Task*>::const_iterator task_it;
			for (task_it = other.tasks_.begin(); task_it!=other.tasks_.end(); task_it++) {
				Task *my_task = new Task();
				*my_task = *(task_it->second);
				tasks_.insert( std::pair<TaskID, Task*>(task_it->first, my_task) );
			}

			// Copy list of higher, lower, equal priority tasks
			higher_prio_tasks_ = other.higher_prio_tasks_;
			lower_prio_tasks_ = other.lower_prio_tasks_;
			equal_prio_tasks_ = other.equal_prio_tasks_;
		}
		
		return *this;
	};
	
	~TaskSet() {
		for (map<TaskID, Task*>::iterator it=tasks_.begin(); it!=tasks_.end(); it++) {
			delete it->second;
		}
	}


 public:	
	map<TaskID, Task*> tasks_; // map task id to its data structure
	map<TaskID, vector<TaskID> > higher_prio_tasks_;
	map<TaskID, vector<TaskID> > lower_prio_tasks_;
	map<TaskID, vector<TaskID> > equal_prio_tasks_;
};

// Function prototypes for generating task set
Task* create_task(int m);
TaskSet* create_taskset(unsigned int m, unsigned int resource_num, unsigned int n_max, CriticalDuration cslen_type);
TaskSet* create_taskset_and_dump(unsigned int m, unsigned int resource_num, unsigned int n_max, 
								 CriticalDuration cslen_type, string parent_folder, bool& is_beginning);

bool init_iteration(TaskSet *taskset, unsigned int m);

#endif // __TASKSET_GENERATOR_H__
