#ifndef __BLOCKING_ANALYSIS_H__
#define __BLOCKING_ANALYSIS_H__

#include "TaskSetGenerator.h"

unsigned int alloc_proc(double C, double L, double D, double Bn, double B1);
double response_time(double C, double L, unsigned int n, double Bn, double B1);
void init_iteration(TaskSet* taskset);
void task_analysis(Task* task, TaskSet* tset);
void blocking_analysis(TaskSet* tset);

#endif
