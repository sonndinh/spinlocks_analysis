#ifndef __TASKSET_GENERATOR_H__
#define __TASKSET_GENERATOR_H__

#include <map>
#include <vector>

using namespace std;

/* Period measured in milliseconds */
#define MIN_TI 10
#define MAX_TI 1000
#define MEAN_TI 500
#define STDEV_TI 200

/* Parameters for short critical section (microsecond) */
#define MIN_SHORT_CSLEN 1
#define MAX_SHORT_CSLEN 15
#define MEAN_SHORT_CSLEN 8
#define STDEV_SHORT_CSLEN 4

/* Parameters for long critical section (microsecond) */
#define MIN_LONG_CSLEN 1
#define MAX_LONG_CSLEN 100
#define MEAN_LONG_CSLEN 50
#define STDEV_LONG_CSLEN 30

/* Moved the definition of Resource sharing factor 
 * to .cpp file to conform to One Definition Rule
 */
//double RSF[] = {0.2, 0.3, 0.5, 0.75};

/* Type of critical section (short or long) */
typedef enum {Short, Long} CriticalDuration;

typedef unsigned int ResourceID;
typedef unsigned int TaskID;

/* Structure for storing resource data under a particular task */
typedef struct Resource {
	ResourceID resourceID;
	double CSLength; // critical section length, in microsecond
	int requestNum;  // number of requests from an enclosed task
} Resource;

/* Implicit deadline parallel task */
typedef struct Task { 
	TaskID taskID;
	double T; // period, in microsecond
	double L; // critical path length, in microsecond
	double U; // utiliation
	double C; // WCET, calculated from u and T
	
	map<ResourceID, Resource*> myResources; // map resouce id to its data

	/* Stored data for fixed-point iteration */
	double B1; // worst-case blocking time on 1 processor, in microsecond
	unsigned int procNum; // number of processors allocated
	double R; // response time, in microsecond
	
	/* Other tasks access to the same resources */
	map<ResourceID, vector<TaskID>*> interferences;

	/* Priority of the task in case of priority lock*/
	unsigned int priority;	
} Task;

typedef struct TaskSet {
	map<TaskID, Task*> tasks; // map task id to its data
} TaskSet;

/* Function prototypes for generating task set */
Task* create_task(int m);
TaskSet* create_taskset(unsigned int m, unsigned int resourceNum, unsigned int Nmax, 
						CriticalDuration cslen_type);

#endif // __TASKSET_GENERATOR_H__
