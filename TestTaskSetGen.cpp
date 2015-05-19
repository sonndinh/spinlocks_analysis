#include <cstdlib>
#include <cstring>
#include <iostream>
#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
	
	if (argc < 4) {
		cout << "Usage: ./prog ProcNum ResourceNum Nmax cs_type" << endl;
		exit(1);
	}

	// Default critical section lengths is short
	CriticalDuration cslen_type = SHORT;
	const int TASKSET_NUM = 1;

	// Number of processors in the system, e.g 16, 32
	const int PROCNUM = atoi(argv[1]);

	// Number of shared resource in the task set, e.g 3, 5
	const int RESOURCE_NUM = atoi(argv[2]);

	// Maximum number of requests to a resource per job, e.g 4
	const int N_MAX = atoi(argv[3]);

	// Set the critical section type
	if ( !strcmp(argv[4], "short") || !strcmp(argv[4], "Short") ) {
		cslen_type = SHORT;
	} else if ( !strcmp(argv[4], "mod") || !strcmp(argv[4], "Mod") ) {
		cslen_type = MODERATE;
	} else if ( !strcmp(argv[4], "long") || !strcmp(argv[4], "Long") ) {
		cslen_type = LONG;
	}

	// Count number of tasksets deemed to be successfully generated
	int success_count = 0;
	// Count number of tasksets deemed to be schedulable after analyzed
	unsigned int count_schedulable = 0;

	/* Generate a bunch of task sets */
	for (int i=0; i<TASKSET_NUM; i++) {
		TaskSet* taskset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		bool ret = init_iteration(taskset, PROCNUM);
		if (ret == true)
			success_count++;

		if ( is_schedulable(taskset, PROCNUM) ) {
			cout << "Taskset is schedulable" << endl;
			count_schedulable++;
		} else {
			cout << "Taskset is unschedulable" << endl;
		}
		/*
		if( blocking_analysis(taskset, PROCNUM) ) {
			cout << "Taskset is schedulable" << endl;
			count_schedulable++;
		} else {
			cout << "Taskset is unschedulable" << endl;
		}
		*/
		
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
		
		map<TaskID, Task*> &tset = taskset->tasks;
		for (map<TaskID, Task*>::iterator it=tset.begin(); it!=tset.end(); it++) {
			free(it->second);
		}
		free(taskset);
	}

	cout << "Percent of schedulable tasksets: " << (double) count_schedulable*100/TASKSET_NUM << "%" << endl;
	cout << "Percent of success after initiating: " << (double)success_count*100/TASKSET_NUM << "%" << endl;
	return 0;
}
