#ifndef __BLOCKING_ANALYSIS_H__
#define __BLOCKING_ANALYSIS_H__

#include <ilcplex/ilocplex.h>
#include <scip/scip.h>
#include "taskset_generator.h"

// Structure to store information of requests 
// from tau_x to l_q which intefere with tau_i
// - requestNum: #requests of tau_x to l_q while tau_i pending
// - csLen: CS length of request from tau_x to l_q (microsecond)
typedef struct {
	unsigned int proc_num;
	unsigned int request_num;
	double cs_len;
} CSData;


class BlockingAnalysis {
 private:
	BlockingAnalysis() {}

 public:
	static BlockingAnalysis& get_instance() {
		static BlockingAnalysis instance;
		return instance;
	}
	
 public:
	// Compute the number of processors required
	static unsigned int alloc_proc(double wcet, double span, double deadline, double whole_blk, double single_blk);
	static double response_time(double wcet, double span, unsigned int proc_num, double whole_blk, double single_blk);

 private:
	// Analyze the upperbound blocking of a task using optimization method
	void task_analysis(Task *task, TaskSet *tset, unsigned int m, SCIP_Real *max_blocking, SpinlockType lock_type);
	unsigned int njobs(Task *task, double t);

	// Generic constraints
	void add_generic_contraint_sum_of_y_for_each_request(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
														 IloRangeArray &cons, Task *task, IloEnv &env);
	void add_generic_constraint_rule_out_mismatched_x_and_y(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
															map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
															IloRangeArray &cons, Task *task, IloEnv &env);
	
	// FIFO-ordered constraints
	void add_fifo_constraint_other_tasks(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
										 map<ResourceID, map<TaskID, CSData> > &x_data, 
										 map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
										 IloRangeArray &cons, Task *task, IloEnv &env);
	void add_fifo_constraint_myself(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
									map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
									IloRangeArray &cons, Task *task, IloEnv &env);
	
	// Priority-ordered common constraints
	void add_prio_constraint_at_most_one_lp_request(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
													map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
													map<ResourceID, map<TaskID, CSData> > &x_data, 
													IloRangeArray &cons, Task *task, TaskSet *taskset, IloEnv &env);
	
	// Priority-orderd with unordered tiebreak constraints
	void add_prio_constraint_unordered_tiebreak(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
												map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
												map<ResourceID, map<TaskID, CSData> > &x_data, 
												IloRangeArray &cons, Task *task, map<ResourceID, double> &dpr, TaskSet *taskset, IloEnv &env);
	void add_prio_constraint_unordered_tiebreak_self(map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
													 map<ResourceID, vector<IloNumVarArray> > &my_x_vars,
													 IloRangeArray &cons, Task *task, IloEnv &env);
	
	// Priority-ordered with FIFO tiebreak constraints
	void add_prio_constraint_fifo_tiebreak_hp(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
											  map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
											  map<ResourceID, map<TaskID, CSData> > &x_data,
											  IloRangeArray &cons, Task *task, map<ResourceID, double> &dpr, TaskSet *taskset, IloEnv &env);
	void add_prio_constraint_fifo_tiebreak_ep(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
											  map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
											  map<ResourceID, map<TaskID, CSData> > &x_data,
											  IloRangeArray &cons, Task *task, TaskSet *taskset, IloEnv &env);
	void add_prio_constraint_fifo_tiebreak_self(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
												map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
												IloRangeArray &cons, Task *task, IloEnv &env);
	
	// Objective function
	void add_objective_function(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
								map<ResourceID, map<TaskID, CSData> > &x_data, 
								map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
								map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
								Task *task, IloModel &model, IloEnv &env);
	
	void call_optimizer(Task *task, TaskSet* taskset, map<ResourceID, map<TaskID, CSData> > &x_data, 
						map<ResourceID, double> &dpr, SCIP_Real *max_blocking, SpinlockType lock_type);
	
	double delay_per_request(Task *task, TaskSet *taskset, ResourceID rid, SpinlockType lock_type);
	SCIP_RETCODE optimize(string fname, SCIP_Real *obj_value);

	// Return upper bound blocking of a task with FIFO locks
	double max_blocking_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data);

	void sort_by_cslen(Task *task, TaskSet *taskset, map<TaskID, CSData> &interfere, vector<CSData> &results);

	// Return upper bound blocking of a task with priority locks
	double max_blocking_prio(Task *task, TaskSet *taskset, map<ResourceID, map<TaskID, CSData> > &x_data,
							 map<ResourceID, double> &dpr, SpinlockType lock_type);

	void blocking_bound_single_proc(Task *task, TaskSet *taskset, unsigned int m, double *blocking, SpinlockType lock_type);

	// Upper bound blocking of a whole job with FIFO locks
	double max_blocking_whole_job_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data);
	void blocking_bound_whole_job(Task *task, TaskSet *taskset, unsigned int m, double *blocking, SpinlockType lock_type);

 public:
	// Schedulability test method
	bool is_schedulable(TaskSet *tset, unsigned int m, SpinlockType lock_type);
};


#endif
