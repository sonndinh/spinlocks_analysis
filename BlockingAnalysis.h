#ifndef __BLOCKING_ANALYSIS_H__
#define __BLOCKING_ANALYSIS_H__

#include <scip/scip.h>
#include "TaskSetGenerator.h"

unsigned int alloc_proc(double C, double L, double D, double Bn, double B1);
double response_time(double C, double L, unsigned int n, double Bn, double B1);
bool init_iteration(TaskSet *taskset, unsigned int m);
void task_analysis(Task *task, TaskSet *tset, unsigned int m, SCIP_Real *max_blocking);
bool blocking_analysis(TaskSet *tset, unsigned int m);
bool is_schedulable(TaskSet *tset, unsigned int m);

/* Helper function prototypes */
unsigned int njobs(Task *task, double t);

#endif
