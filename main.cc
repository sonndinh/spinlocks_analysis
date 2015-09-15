#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include "taskset_generator.h"
#include "blocking_analysis.h"

using namespace std;


// Global variables accessed by other files
int g_proc_num;
int g_resource_num;
int g_max_request_num;
string g_cs_type;
int g_taskset_num;

// Only use for debugging: return true if 2 tasksets are equal
// otherwise, return false
bool compare_tasksets(const TaskSet &first, const TaskSet &second) {
	assert(first.tasks_.size() == second.tasks_.size());
	for (map<TaskID, Task*>::const_iterator it=first.tasks_.begin(); it!=first.tasks_.end(); it++) {
		Task* first_task = it->second;
		Task* second_task = second.tasks_.at(it->first);
		assert(second_task != NULL);
		assert(first_task->period_ == second_task->period_);
		assert(first_task->span_ == second_task->span_);
		assert(first_task->utilization_ == second_task->utilization_);
		assert(first_task->wcet_ == second_task->wcet_);
		assert(first_task->priority_ == second_task->priority_);
		assert(first_task->blocking_on_single_ == second_task->blocking_on_single_);
		assert(first_task->proc_num_ == second_task->proc_num_);
		assert(first_task->response_time_ == second_task->response_time_);
		
		assert(first_task->my_resources_.size() == second_task->my_resources_.size());
		for(map<ResourceID, Resource*>::iterator iter=first_task->my_resources_.begin(); 
			iter != first_task->my_resources_.end(); iter++) {
			Resource* first_res = iter->second;
			Resource* second_res = second_task->my_resources_.at(iter->first);
			assert(second_res != NULL);
			assert(first_res->resource_id == second_res->resource_id);
			assert(first_res->critical_section_len == second_res->critical_section_len);
			assert(first_res->request_num == second_res->request_num);
		}

		assert(first_task->interferences_.size() == second_task->interferences_.size());
		for (map<ResourceID, vector<TaskID>*>::iterator iter=first_task->interferences_.begin();
			 iter != first_task->interferences_.end(); iter++) {
			vector<TaskID>* first_vec = iter->second;
			vector<TaskID>* second_vec = second_task->interferences_.at(iter->first);
			assert(second_vec != NULL);
			assert(first_vec->size() == second_vec->size());
			for (int i=0; i<first_vec->size(); i++) {
				assert(first_vec->at(i) == second_vec->at(i));
			}
		}
	}

	
	for (map<TaskID, vector<TaskID> >::const_iterator it=first.higher_prio_tasks_.begin();
		 it != first.higher_prio_tasks_.end(); it++) {
		const vector<TaskID> &first_vec = it->second;
		const vector<TaskID> &second_vec = second.higher_prio_tasks_.at(it->first);
		assert(first_vec.size() == second_vec.size());
		for (int i=0; i<first_vec.size(); i++) {
			assert(first_vec[i] == second_vec[i]);
		}
	}

	for (map<TaskID, vector<TaskID> >::const_iterator it=first.lower_prio_tasks_.begin();
		 it != first.lower_prio_tasks_.end(); it++) {
		const vector<TaskID> &first_vec = it->second;
		const vector<TaskID> &second_vec = second.lower_prio_tasks_.at(it->first);
		assert(first_vec.size() == second_vec.size());
		for (int i=0; i<first_vec.size(); i++) {
			assert(first_vec[i] == second_vec[i]);
		}
	}

	for (map<TaskID, vector<TaskID> >::const_iterator it=first.equal_prio_tasks_.begin();
		 it != first.equal_prio_tasks_.end(); it++) {
		const vector<TaskID> &first_vec = it->second;
		const vector<TaskID> &second_vec = second.equal_prio_tasks_.at(it->first);
		assert(first_vec.size() == second_vec.size());
		for (int i=0; i<first_vec.size(); i++) {
			assert(first_vec[i] == second_vec[i]);
		}
	}
	
	return true;
}

int main(int argc, char* argv[]) {
	
	if (argc < 6) {
		cout << "Usage: " << argv[0] 
			 << " {processor_number} {resource_number} {max_critical_sections_number} {critical_section_type} {repetition}" << endl;
		exit(1);
	}

	// Default critical section lengths is short
	CriticalDuration cslen_type;

	// Number of processors in the system, e.g 16, 32
	const int kProcNum = atoi(argv[1]);
	g_proc_num = kProcNum;

	// Number of shared resource in the task set, e.g 3, 5
	const int kResourceNum = atoi(argv[2]);
	g_resource_num = kResourceNum;

	// Maximum number of requests to a resource per job, e.g 4
	const int kNMax = atoi(argv[3]);
	g_max_request_num = kNMax;

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

	const int kTasksetNum = atoi(argv[5]);
	g_taskset_num = kTasksetNum;

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

	// Beginning of the run
	bool is_beginning = true;
	BlockingAnalysis& analyzer = BlockingAnalysis::get_instance();
	// Generate a bunch of task sets
	for (int i=0; i<kTasksetNum; i++) {
		bool ret;
		//TaskSet* taskset1 = create_taskset(kProcNum, kResourceNum, kNMax, cslen_type);
		string parent_folder = "/export/austen/home/son/codes/spinlocks_analysis/tasksets";
		TaskSet* taskset1 = create_taskset_and_dump(kProcNum, kResourceNum, kNMax, cslen_type, parent_folder, is_beginning);
		ret = init_iteration(taskset1, kProcNum);
		if (ret == true)
			success_count++;

		if ( analyzer.is_schedulable(taskset1, kProcNum, FIFO) ) {
			count_schedulable_fifo++;
		}

		// Generate another task set and analyze with priority-ordered unordered tiebreak
		//		TaskSet *another_tset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		//		ret = init_iteration(another_tset, PROCNUM);
		//		TaskSet *taskset2 = new TaskSet();
		//		*taskset2 = *taskset1;

		// Test the copy assignment operator of TaskSet
		//		assert(compare_tasksets(*taskset1, *taskset2));
		
		//		cout << "--- PRIORITY LOCKS: " << endl;
		//		if ( analyzer.is_schedulable(taskset2, kProcNum, PRIO_UNORDERED) ) {
			//			cout << "Taskset is schedulable with Priority-ordered spin locks" << endl;
		//			count_schedulable_prio++;
		//		} else {
			//			cout << "Taskset is unschedulable with Priority-ordered spin locks" << endl;
		//		}

		// Generate another task set and analyze with priority-ordered fifo tiebreak
		//		TaskSet *another_tset2 = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		//		ret = init_iteration(another_tset2, PROCNUM);
		TaskSet *taskset3 = new TaskSet();
		*taskset3 = *taskset1;

		assert(compare_tasksets(*taskset1, *taskset3));
		
		if ( analyzer.is_schedulable(taskset3, kProcNum, PRIO_FIFO) ) {
			count_schedulable_prio_fifo++;
		}

		delete taskset1;
		//delete taskset2;
		delete taskset3;
	}

	/*
	fifo_ofile << (double) count_schedulable_fifo/TASKSET_NUM;
	prio_unordered_ofile << (double) count_schedulable_prio/TASKSET_NUM;
	prio_fifo_ofile << (double) count_schedulable_prio_fifo/TASKSET_NUM;
	fifo_ofile.close();
	prio_unordered_ofile.close();
	prio_fifo_ofile.close();	
	*/
	
	cout << "Percent of schedulable tasksets with FIFO locks: " << (double) count_schedulable_fifo*100/kTasksetNum << "%" << endl;
	//	cout << "Percent of schedulable tasksets with Priority locks: " << (double) count_schedulable_prio*100/kTasksetNum << "%" << endl;
	cout << "Percent of schedulable tasksets with Priority locks & FIFO tiebreak: " << (double) count_schedulable_prio_fifo*100/kTasksetNum << "%" << endl;

	return 0;
}
