#include "TaskSetGenerator.h"
#include "BlockingAnalysis.h"

int main(int argc, char** argv) {
	CriticalDuration cslen_type = Short;
	TaskSet* tset = create_taskset(16, 8, 5, cslen_type);
	init_iteration(tset);
	return 0;
}
