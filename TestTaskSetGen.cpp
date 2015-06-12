#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
	
	if (argc != 6) {
		cout << "Usage: ./prog ProcNum ResourceNum Nmax cs_type repetition" << endl;
		exit(1);
	}

	// Default critical section lengths is short
	CriticalDuration cslen_type;
	//	const int TASKSET_NUM = 10;

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
	} else if (argv[4] == NULL) {
		cslen_type = SHORT;
	}

	const int TASKSET_NUM = atoi(argv[5]);

	// Count number of tasksets deemed to be successfully generated
	int success_count = 0;
	// Count number of tasksets deemed to be schedulable after analyzed
	unsigned int count_schedulable_fifo = 0;
	unsigned int count_schedulable_prio = 0;
	unsigned int count_schedulable_prio_fifo = 0;

	// Output to a file
	ofstream outfile;
	stringstream ss;
	//	ss << "results/exp_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_" << argv[4] << "_" << argv[5] << ".dat";
	//	outfile.open(ss.str().c_str());
	
	// Generate a bunch of task sets
	for (int i=0; i<TASKSET_NUM; i++) {
		TaskSet* taskset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		bool ret = init_iteration(taskset, PROCNUM);
		if (ret == true)
			success_count++;

		if ( is_schedulable2(taskset, PROCNUM, FIFO) ) {
			cout << "Taskset is schedulable with FIFO spin locks" << endl;
			count_schedulable_fifo++;
		} else {
			cout << "Taskset is unschedulable with FIFO spin locks" << endl;
		}

		// Generate another task set and analyze with priority-ordered unordered tiebreak
		TaskSet *another_tset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		ret = init_iteration(another_tset, PROCNUM);
		
		if ( is_schedulable2(another_tset, PROCNUM, PRIO_UNORDERED) ) {
			cout << "Taskset is schedulable with Priority-ordered spin locks" << endl;
			count_schedulable_prio++;
		} else {
			cout << "Taskset is unschedulable with Priority-ordered spin locks" << endl;
		}

		// Generate another task set and analyze with priority-ordered fifo tiebreak
		TaskSet *another_tset2 = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		ret = init_iteration(another_tset2, PROCNUM);
		
		if (is_schedulable2(another_tset2, PROCNUM, PRIO_FIFO) ) {
			cout << "Taskset is schedulable with Prioirty-ordered spin locks & FIFO tiebreak" << endl;
			count_schedulable_prio_fifo++;
		} else {
			cout << "Taskset is unschedulable with Priority-ordered spin locks & FIFO tiebreak" << endl;
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

	//	outfile << (double) count_schedulable*100/TASKSET_NUM;
	//	outfile.close();
	
	cout << "Percent of schedulable tasksets with FIFO locks: " << (double) count_schedulable_fifo*100/TASKSET_NUM << "%" << endl;
	cout << "Percent of schedulable tasksets with Priority locks: " << (double) count_schedulable_prio*100/TASKSET_NUM << "%" << endl;
	cout << "Percent of schedulable tasksets with Priority locks & FIFO tiebreak: " << (double) count_schedulable_prio_fifo*100/TASKSET_NUM << "%" << endl;
	//	cout << "Percent of success after initiating: " << (double)success_count*100/TASKSET_NUM << "%" << endl;
	return 0;
}
