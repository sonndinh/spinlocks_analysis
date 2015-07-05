#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

using namespace std;


// Global variables accessed by other files
int g_proc_num;
int g_resource_num;
int g_max_request_num;
string g_cs_type;
int g_taskset_num;


extern int large_critical_path_length_fail_num;
extern int different_num_opt_numerical;
extern int equal_num_opt_numerical;
extern int opt_larger_numerical_num;
extern int numerical_larger_opt_num;
int main(int argc, char** argv) {
	
	if (argc < 6) {
		cout << "Usage: ./prog ProcNum ResourceNum Nmax cs_type repetition" << endl;
		exit(1);
	}

	// Default critical section lengths is short
	CriticalDuration cslen_type;
	//	const int TASKSET_NUM = 10;

	// Number of processors in the system, e.g 16, 32
	const int PROCNUM = atoi(argv[1]);
	g_proc_num = PROCNUM;

	// Number of shared resource in the task set, e.g 3, 5
	const int RESOURCE_NUM = atoi(argv[2]);
	g_resource_num = RESOURCE_NUM;

	// Maximum number of requests to a resource per job, e.g 4
	const int N_MAX = atoi(argv[3]);
	g_max_request_num = N_MAX;

	// Set the critical section type
	if ( !strcmp(argv[4], "short") || !strcmp(argv[4], "Short") ) {
		cslen_type = SHORT;
		g_cs_type = "short";
	} else if ( !strcmp(argv[4], "mod") || !strcmp(argv[4], "Mod") ) {
		cslen_type = MODERATE;
		g_cs_type = "mod";
	} else if ( !strcmp(argv[4], "long") || !strcmp(argv[4], "Long") ) {
		cslen_type = LONG;
		g_cs_type = "long";
	} else {
		cout << "Wrong critical section type!" << endl;
		exit(1);
	}

	const int TASKSET_NUM = atoi(argv[5]);
	g_taskset_num = TASKSET_NUM;

	// Count number of tasksets deemed to be successfully generated
	int success_count = 0;
	// Count number of tasksets deemed to be schedulable after analyzed
	unsigned int count_schedulable_fifo = 0;
	unsigned int count_schedulable_prio = 0;
	unsigned int count_schedulable_prio_fifo = 0;

	// Output to a file
	ofstream fifo_ofile;
	ofstream prio_unordered_ofile;
	ofstream prio_fifo_ofile;
	stringstream ss;
	stringstream header;
	header << g_proc_num << "_" << g_resource_num << "_" << g_max_request_num << "_" << g_cs_type << "_" << g_taskset_num;

	//	ss << "results/exp_" << argv[1] << "_" << argv[2] << "_" << argv[3] << "_" << argv[4] << "_" << argv[5] << ".dat";

	/*
	// Open output file for FIFO case
	ss << "results/" << header.str() << "_fifo" << ".dat";
	fifo_ofile.open(ss.str().c_str());
	ss.str(string());

	// Open output file for Priority with unordered tiebreak
	ss << "results/" << header.str() << "_prio" << ".dat";
	prio_unordered_ofile.open(ss.str().c_str());
	ss.str(string());

	// Open output file for Priority with FIFO tiebreak
	ss << "results/" << header.str() << "_prio_fifo" << ".dat";
	prio_fifo_ofile.open(ss.str().c_str());
	*/
	
	// Generate a bunch of task sets
	for (int i=0; i<TASKSET_NUM; i++) {
		bool ret;
		/*
		TaskSet* taskset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		ret = init_iteration(taskset, PROCNUM);
		if (ret == true)
			success_count++;

		if ( is_schedulable2(taskset, PROCNUM, FIFO) ) {
			cout << "Taskset is schedulable with FIFO spin locks" << endl;
			count_schedulable_fifo++;
		} else {
			cout << "Taskset is unschedulable with FIFO spin locks" << endl;
		}
		*/

		// Generate another task set and analyze with priority-ordered unordered tiebreak
		TaskSet *another_tset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		ret = init_iteration(another_tset, PROCNUM);
		
		if ( is_schedulable2(another_tset, PROCNUM, PRIO_UNORDERED) ) {
			cout << "Taskset is schedulable with Priority-ordered spin locks" << endl;
			count_schedulable_prio++;
		} else {
			cout << "Taskset is unschedulable with Priority-ordered spin locks" << endl;
		}

		/*
		// Generate another task set and analyze with priority-ordered fifo tiebreak
		TaskSet *another_tset2 = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		ret = init_iteration(another_tset2, PROCNUM);
		
		if (is_schedulable2(another_tset2, PROCNUM, PRIO_FIFO) ) {
			cout << "Taskset is schedulable with Prioirty-ordered spin locks & FIFO tiebreak" << endl;
			count_schedulable_prio_fifo++;
		} else {
			cout << "Taskset is unschedulable with Priority-ordered spin locks & FIFO tiebreak" << endl;
		}
		*/

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
		
		/*
		map<TaskID, Task*> &tset = taskset->tasks;
		for (map<TaskID, Task*>::iterator it=tset.begin(); it!=tset.end(); it++) {
			free(it->second);
		}
		free(taskset);
		*/
	}

	/*
	fifo_ofile << (double) count_schedulable_fifo/TASKSET_NUM;
	prio_unordered_ofile << (double) count_schedulable_prio/TASKSET_NUM;
	prio_fifo_ofile << (double) count_schedulable_prio_fifo/TASKSET_NUM;
	fifo_ofile.close();
	prio_unordered_ofile.close();
	prio_fifo_ofile.close();	
	*/
	
	//	cout << "Percent of schedulable tasksets with FIFO locks: " << (double) count_schedulable_fifo*100/TASKSET_NUM << "%" << endl;
	cout << "Percent of schedulable tasksets with Priority locks: " << (double) count_schedulable_prio*100/TASKSET_NUM << "%" << endl;
	//cout << "Percent of schedulable tasksets with Priority locks & FIFO tiebreak: " << (double) count_schedulable_prio_fifo*100/TASKSET_NUM << "%" << endl;

	//	cout << "The number of failures caused by large inflated critical path length: " << large_critical_path_length_fail_num << endl;
	cout << "Number of times blocking by optimization and numerical method are different: " << different_num_opt_numerical << endl;
	cout << "Number of times blocking by optimization > by numerical: " << opt_larger_numerical_num << endl;
	cout << "Number of times blocking by numerical > by optimization: " << numerical_larger_opt_num << endl;
	cout << "Number of times blocking by optimization and numerical method are same: " << equal_num_opt_numerical << endl;
	//	cout << "Percent of success after initiating: " << (double)success_count*100/TASKSET_NUM << "%" << endl;
	return 0;
}
