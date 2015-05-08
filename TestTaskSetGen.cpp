#include <cstdlib>
#include <iostream>
#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
	
	if (argc != 4) {
		cout << "Usage: ./prog ProcNum ResourceNum Nmax" << endl;
		exit(1);
	}
	
	CriticalDuration cslen_type = Short;
	const int TASKSET_NUM = 10;
	// Number of processors in the system, e.g 16, 32
	const int PROCNUM = atoi(argv[1]);
	// Number of shared resource in the task set, e.g 3, 5
	const int RESOURCE_NUM = atoi(argv[2]);
	// Maximum number of requests to a resource per job, e.g 4
	const int N_MAX = atoi(argv[3]);
	int success_count = 0;

	unsigned int count_schedulable = 0;
	/* Generate a bunch of task sets */
	for (int i=0; i<TASKSET_NUM; i++) {
		TaskSet* taskset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		bool ret = init_iteration(taskset, PROCNUM);
		if (ret == true)
			success_count++;
		
		if( blocking_analysis(taskset, PROCNUM) ) {
			cout << "Taskset is schedulable" << endl;
			count_schedulable++;
		} else {
			cout << "Taskset is unschedulable" << endl;
		}
		
		/*
		map<TaskID, Task*> &tset = taskset->tasks;
		map<TaskID, Task*>::iterator it = tset.begin();
		for (; it!=tset.end(); it++) {
			//			if (it->first == 2) {
				// test task 2 only
				task_analysis(it->second, taskset, PROCNUM);
				//	break;
				//			}
		}
		*/
	}

	cout << "Percent of schedulable tasksets: " << (double) count_schedulable*100/TASKSET_NUM << "%" << endl;
	//	cout << "Percent of success after initiating: " << (double)success_count*100/TASKSET_NUM << "%" << endl;
	return 0;
}
