#include <iostream>
#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

using namespace std;

int main(int argc, char** argv) {
	CriticalDuration cslen_type = Short;
	const int TASKSET_NUM = 1;
	const int PROCNUM = 16;
	const int RESOURCE_NUM = 3;
	const int N_MAX = 4;
	int success_count = 0;

	/* Generate a bunch of task sets */
	for (int i=0; i<TASKSET_NUM; i++) {
		TaskSet* taskset = create_taskset(PROCNUM, RESOURCE_NUM, N_MAX, cslen_type);
		bool ret = init_iteration(taskset, PROCNUM);
		if (ret == true)
			success_count++;
		
		map<TaskID, Task*> &tset = taskset->tasks;
		map<TaskID, Task*>::iterator it = tset.begin();
		for (; it!=tset.end(); it++) {
			task_analysis(it->second, taskset, PROCNUM);
		}
	}

	cout << "Percent of success after initiating: " << (double)success_count*100/TASKSET_NUM << "%" << endl;
	return 0;
}
