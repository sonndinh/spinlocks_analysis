#include <cmath>
#include <iostream>
#include <algorithm>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include "blocking_analysis.h"

ILOSTLBEGIN


extern int g_proc_num;
extern int g_resource_num;
extern int g_max_request_num;
extern string g_cs_type;
extern int g_taskset_num;

// Update number of processors allocated
unsigned int BlockingAnalysis::alloc_proc(double wcet, double span, double deadline, double whole_blk, double single_blk) {
	if (deadline < span + single_blk)
		cout << "Denominator of equation for n_i is negative !!!" << endl;
	else if (deadline == span + single_blk)
		cout << "Denominator of equation for n_i is zero !!!" << endl;

	return ceil(((wcet+whole_blk) - (span+single_blk))/(deadline - (span+single_blk)));
}

// Update bound of response time using blocking time, #processors
double BlockingAnalysis::response_time(double wcet, double span, unsigned int proc_num, double whole_blk, double single_blk) {
	double response_time = (wcet+whole_blk-span-single_blk)/proc_num + span + single_blk;
	return ceil(response_time); // in microsecond
}

// Update task's blocking time, #processors, and response time 
// based on results from the previous iteration.
void BlockingAnalysis::task_analysis(Task* task, TaskSet* taskset, unsigned int m, SCIP_Real *max_blocking, SpinlockType lock_type) {
	TaskID my_id = task->task_id_;
	TaskID r_i = task->response_time_;
	map<ResourceID, vector<TaskID>*> &interferences = task->interferences_;
	map<ResourceID, vector<TaskID>*>::iterator rit = interferences.begin();
	map<TaskID, Task*> &tset = taskset->tasks_;

	// Gather information from requests of tau_x which 
	// interfere with my requests
	map<ResourceID, map<TaskID, CSData> > x_data;
	for (; rit != interferences.end(); rit++) {
		ResourceID res_id = rit->first;
		vector<TaskID>* vec = rit->second;
		for (int i=0; i<vec->size(); i++) {
			TaskID task_id = vec->at(i);
			// Abort if this is me, but this is impossible
			if (task_id == my_id)
				continue;
			Task* tau_x = tset[task_id];
			unsigned int proc_num = tau_x->proc_num_;
			map<ResourceID, Resource*> &tau_x_resources = tau_x->my_resources_;
			unsigned int n_x_q = tau_x_resources[res_id]->request_num;
			double cs_len = tau_x_resources[res_id]->critical_section_len;
			unsigned int request_num = njobs(tau_x, r_i)*n_x_q;
			map<TaskID, CSData> &inner_map = x_data[res_id];
			CSData data = {proc_num, request_num, cs_len};
			inner_map.insert(std::pair<TaskID, CSData> (task_id, data));
		}
	}

	/*
	// Print debug information
	//#define _ITER_DEBUG_
#ifdef _ITER_DEBUG_
	cout << "FOR TASK " << myId << endl;
	map<ResourceID, map<TaskID, CSData> >::iterator it = x_data.begin();
	for (; it != x_data.end(); it++) {
		ResourceID rId = it->first;
		cout << "Resource ID: " << rId << endl;
		cout << "Interfere (TaskID, RequestNum, CSLen): ";
		map<TaskID, CSData> &inner_map = it->second;
		map<TaskID, CSData>::iterator inner_it = inner_map.begin();
		for (; inner_it != inner_map.end(); inner_it++) {
			cout << "(" << inner_it->first << "," <<
				inner_it->second.requestNum << "," <<
				inner_it->second.csLen << ")  ";
		}
		cout << endl;
	}
	cout << endl;
#endif  // _ITER_DEBUG_
	*/

	// Calculate delay-per-request if lock type is priority-ordered
	map<ResourceID, double> dpr;
	if (lock_type == PRIO_UNORDERED || lock_type == PRIO_FIFO) {
		map<ResourceID, Resource*> &my_resources = task->my_resources_;
		for (map<ResourceID, Resource*>::iterator it = my_resources.begin(); it != my_resources.end(); it++) {
			ResourceID res_id = it->first;
			double delay = delay_per_request(task, taskset, res_id, lock_type);
			dpr[res_id] = delay;
		}
	}
	
	// Call the optimizer to solve the maximum blocking time
	SCIP_Real max_blk;
	call_optimizer(task, taskset, x_data, dpr, &max_blk, lock_type);
	*max_blocking = max_blk;
}


void BlockingAnalysis::call_optimizer(Task *task, TaskSet *taskset, map<ResourceID, map<TaskID, CSData> > &x_data, 
									  map<ResourceID, double> &dpr, SCIP_Real *blocking, SpinlockType lock_type) {
	// Solve the maximization of blocking
	IloEnv env;
	try {
		IloModel model(env);
		IloRangeArray cons(env);
		// Variables table of requests from other tasks
		map<ResourceID, map<TaskID, IloNumVarArray> > x_vars;
		map<ResourceID, map<TaskID, CSData> >::iterator iter = x_data.begin();

		for (; iter != x_data.end(); iter++) {
			ResourceID res_id = iter->first;

			// Get an inner map for this resource ID
			map<TaskID, IloNumVarArray> &inner_map = x_vars[res_id];
			map<TaskID, CSData> &task_data = iter->second;
			map<TaskID, CSData>::iterator inner_iter = task_data.begin();
			for (; inner_iter != task_data.end(); inner_iter++) {
				// Each other task has a corresponding IloNumVarArray
				TaskID task_id = inner_iter->first;
				// Get number of interfering requests from this task 
				// to this resource
				unsigned int var_num = inner_iter->second.request_num;
				IloNumVarArray x(env);
				// All x variables are in range [0,1]
				for (int i=0; i<var_num; i++) {
					// Trick: Use a very small number (instead of 0) for lower bound of x vars
					//					x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
					// Use 0 as lower bound of x vars
					x.add(IloNumVar(env, 0, 1.0, ILOFLOAT));
				}
				inner_map.insert(std::pair<TaskID, IloNumVarArray>(task_id, x));
			} // Task
		} // Resource
		
		// Populate data for my (tau_i's) y & x variables 
		// For my_y_vars: each element of a vector 
		// corresponds to a row of matrix Y (n_i * N_i_q)
		// For my_x_vars: each element of a vector 
		// corresponds to a column of matrix X (N_i_q * n_i)
		unsigned int my_proc_num = task->proc_num_;
		map<ResourceID, vector<IloNumVarArray> > my_y_vars;
		map<ResourceID, vector<IloNumVarArray> > my_x_vars;
		map<ResourceID, Resource*> &my_res = task->my_resources_;
		map<ResourceID, Resource*>::iterator my_rit = my_res.begin();
		
		for (; my_rit != my_res.end(); my_rit++) {
			ResourceID r_id = my_rit->first;
			unsigned int my_req_num = my_rit->second->request_num;
			vector<IloNumVarArray> &ys = my_y_vars[r_id];
			vector<IloNumVarArray> &xs = my_x_vars[r_id];
			for (int u=0; u<my_proc_num; u++) {
				IloNumVarArray y(env);
				IloNumVarArray x(env);
				for (int k=0; k<my_req_num; k++) {
					y.add(IloNumVar(env, 0, 1, ILOINT));
					// Trick: use a very small number for lower bound of x vars
					//					x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
					// Use 0 as lower bound of x vars
					x.add(IloNumVar(env, 0, 1.0, ILOFLOAT));
				}
				ys.push_back(y);
				xs.push_back(x);
			}
		}

		// Constraint 3: sum of each y column is <= 1
		add_generic_contraint_sum_of_y_for_each_request(my_y_vars, cons, task, env);

		// Constraint 5: mismatched product of y&x is zero
		add_generic_constraint_rule_out_mismatched_x_and_y(my_y_vars, my_x_vars, cons, task, env);

		if (lock_type == FIFO) {
			// Constraint 8.1: at most 1 request from each other processor of another task block me
			add_fifo_constraint_other_tasks(x_vars, x_data, my_y_vars, cons, task, env);

			// Constraint 8.2: doing the same for other processors in the same task
			add_fifo_constraint_myself(my_y_vars, my_x_vars, cons, task, env);

		} else if (lock_type == PRIO_UNORDERED) {
			// Constraint 9: at most 1 lower locking-priority request can precede each of my requests
			add_prio_constraint_at_most_one_lp_request(x_vars, my_y_vars, x_data, cons, task, taskset, env);
			
			// Constraint 10: bound contribution of equal & higher locking-priority tasks
			add_prio_constraint_unordered_tiebreak(x_vars, my_y_vars, x_data, cons, task, dpr, taskset, env);
			
			// Constraint for blocking caused by other processors of the same task
			add_prio_constraint_unordered_tiebreak_self(my_y_vars, my_x_vars, cons, task, env);
			
		} else if (lock_type == PRIO_FIFO) {
			// Constraint 9: apply to both types of priority-ordered spin locks
			add_prio_constraint_at_most_one_lp_request(x_vars, my_y_vars, x_data, cons, task, taskset, env);

			// Constraint 11: bound contribution of higher locking-priority tasks
			add_prio_constraint_fifo_tiebreak_hp(x_vars, my_y_vars, x_data, cons, task, dpr, taskset, env);

			// Constraint 12: bound contribution of equal locking-priotity tasks
			add_prio_constraint_fifo_tiebreak_ep(x_vars, my_y_vars, x_data, cons, task, taskset, env);
			
			// Constraint 13: bound contribution of the other processors of the same task
			add_prio_constraint_fifo_tiebreak_self(my_y_vars, my_x_vars, cons, task, env);
		}
		
		// Add constraints to Cplex model
		model.add(cons);
		
		// Formulate the objective function & add obj function to Cplex model
		add_objective_function(x_vars, x_data, my_y_vars, my_x_vars, task, model, env);

		IloCplex cplex(model);

		// disable output to stdout
		cplex.setOut(env.getNullStream());

		// create file name
		TaskID task_id = task->task_id_;
		string file_name;
		stringstream ss; // full path
		stringstream header; // inner folder name
		header << g_proc_num << "_" << g_resource_num << "_" << g_max_request_num << "_" << g_cs_type;
		if (lock_type == FIFO) {
			ss << "tmp/" << header.str() << "_fifo_task_" << task_id << ".lp";
		} else if (lock_type == PRIO_UNORDERED) {
			ss << "tmp/" << header.str() << "_prio_task_" << task_id << ".lp";
		} else if (lock_type == PRIO_FIFO) {
			ss << "tmp/" << header.str() << "_prio_fifo_task_" << task_id << ".lp";
		}
		file_name = ss.str();

		// export to .lp file
		cplex.exportModel(file_name.c_str());
		
		SCIP_Real sol_val;
		// SCIP read & solve input from the .lp file
		SCIP_RETCODE retcode = optimize(file_name, &sol_val);
		if (retcode != SCIP_OKAY) {
			SCIPprintError(retcode);
			exit(-1);
		}
		*blocking = sol_val;


	} catch (IloException &e) {
		cerr << "Ilog concert exception: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception: " << endl;
	}

	cout << endl;
	env.end();
}

SCIP_RETCODE BlockingAnalysis::optimize(string fname, SCIP_Real *obj_value) {
	const char* file_name= fname.c_str();
	SCIP *scip = NULL;
	SCIP_CALL( SCIPcreate(&scip) );
	
	// load default plugins
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
	
	// disable output to stdout
	SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );

	// Set time limit to 300 seconds (5 minutes)
	//	SCIP_Real timelimit;
	//	SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
	//	cout << "Default time limit of SCIP: " << timelimit << endl;
	SCIP_CALL( SCIPsetRealParam(scip, "limits/time", 300) );
	
	SCIP_CALL( SCIPreadProb(scip, file_name, NULL) );
	
	SCIP_CALL( SCIPsolve(scip) );
	
	// get best solution
	SCIP_Real obj_val = SCIPgetPrimalbound(scip);
	std::cout << "Objective value: " << obj_val << std::endl;

	*obj_value = obj_val;
	
	SCIP_CALL( SCIPfree(&scip) );
	
	BMScheckEmptyMemory();
	return SCIP_OKAY;
}

// Constraint 3: for each request to a resource, sum 
// of corresponding y variables is at most 1
// Constraint 2 then follows
void BlockingAnalysis::add_generic_contraint_sum_of_y_for_each_request(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
																	   IloRangeArray &cons, Task *task, IloEnv &env) {	
	unsigned int my_proc_num = task->proc_num_;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	map<ResourceID, Resource*> &my_res = task->my_resources_;
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID r_id = y_it->first;			
		vector<IloNumVarArray> &var_arrays = y_it->second;
		unsigned int req_num = my_res[r_id]->request_num;
		for (int k=0; k<req_num; k++) {
			IloExpr expr(env);
			for (int u=0; u<my_proc_num; u++) {
				expr += var_arrays[u][k];
			}
			cons.add(expr <= 1);
			expr.end(); // Important: must free it
		}
	}
}

// Constraint 5: products of mismatched y and x variables 
// are all zero
void BlockingAnalysis::add_generic_constraint_rule_out_mismatched_x_and_y(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
																		  map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
																		  IloRangeArray &cons, Task *task, IloEnv &env) {
	unsigned int my_proc_num = task->proc_num_;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	map<ResourceID, Resource*> &my_res = task->my_resources_;
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID r_id = y_it->first;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[r_id];
		unsigned int req_num = my_res[r_id]->request_num;
		for (int yi=0; yi<my_proc_num; yi++) {
			for (int xi=0; xi<my_proc_num; xi++) {
				if (yi == xi)
					continue;
				
				IloExpr expr(env);
				for (int k=0; k<req_num; k++) {
					expr += ys[yi][k] * xs[xi][k];
				}
				cons.add(expr <= 0);
				expr.end();
			}
		}
		
		// Constraint 6: requests from me cannot block me 
		// (processor under consideration is the first 
		// processor of tau_i)
		IloExpr expr(env);
		for (int k=0; k<req_num; k++) {
			expr += ys[0][k] * xs[0][k];
		}
		cons.add(expr <= 0);
		expr.end();
	}
}


// Constraint 8 (FIFO): at most 1 request from each other
// processor of another task can block a single request from me
void BlockingAnalysis::add_fifo_constraint_other_tasks(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
													   map<ResourceID, map<TaskID, CSData> > &x_data, 
													   map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
													   IloRangeArray &cons, Task *task, IloEnv &env) {
	map<ResourceID, Resource*> &my_res = task->my_resources_;
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator inter_it = x_vars.begin();
	for (; inter_it != x_vars.end(); inter_it++) {
		ResourceID r_id = inter_it->first;
		map<TaskID, IloNumVarArray> &others = inter_it->second;
		map<TaskID, IloNumVarArray>::iterator others_it = others.begin();
		for (; others_it != others.end(); others_it++) {
			TaskID t_id = others_it->first;
			IloNumVarArray &x = others_it->second;
			IloExpr expr(env);
			unsigned int his_request_num = x_data[r_id][t_id].request_num;
			int his_proc_num = x_data[r_id][t_id].proc_num;
			for (int k=0; k<his_request_num; k++) {
				expr += x[k];
			}
			vector<IloNumVarArray> &my_ys = my_y_vars[r_id];
			unsigned int my_request_num = my_res[r_id]->request_num;
			for (int i=0; i<my_request_num; i++) {
				expr -= his_proc_num * my_ys[0][i];
			}
			cons.add(expr <= 0);
			expr.end();
		}
	}
}

// Constraint 8 (FIFO): for interference within myself
void BlockingAnalysis::add_fifo_constraint_myself(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
												  map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
												  IloRangeArray &cons, Task *task, IloEnv &env) {
	unsigned int my_proc_num = task->proc_num_;
	map<ResourceID, Resource*> &my_res = task->my_resources_;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID r_id = y_it->first;
		unsigned int my_request_num = my_res[r_id]->request_num;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[r_id];
		assert(my_proc_num == ys.size());
		IloExpr rhs(env);
		for (int i=0; i<my_request_num; i++) {
			rhs += ys[0][i];
		}
		
		for (int i=1; i<my_proc_num; i++) {
			IloExpr lhs(env);
			for (int j=0; j<my_request_num; j++) {
				lhs += ys[i][j] * xs[i][j];
			}
			cons.add(lhs - rhs <= 0);
			lhs.end();
		}
		rhs.end();
	}
}

// Build objective function for the blocking time
void BlockingAnalysis::add_objective_function(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
											  map<ResourceID, map<TaskID, CSData> > &x_data, 
											  map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
											  map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
											  Task *task, IloModel &model, IloEnv &env) {
	IloExpr obj(env);
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator inter_it = x_vars.begin();
	for (; inter_it != x_vars.end(); inter_it++) {
		ResourceID r_id = inter_it->first;
		map<TaskID, IloNumVarArray> &others = inter_it->second;
		map<TaskID, IloNumVarArray>::iterator others_it = others.begin();
		for (; others_it != others.end(); others_it++) {
			TaskID t_id = others_it->first;
			double cs_len = x_data[r_id][t_id].cs_len;
			unsigned int var_num = x_data[r_id][t_id].request_num;
			for (int i=0; i<var_num; i++) {
				obj += cs_len * others[t_id][i];
			}
		}
	}
	
	unsigned int my_proc_num = task->proc_num_;
	map<ResourceID, Resource*> &my_res = task->my_resources_;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID r_id = y_it->first;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[r_id];
		unsigned int my_req_num = my_res[r_id]->request_num;
		double cs_len = my_res[r_id]->critical_section_len;
		
		for (int u=0; u<my_proc_num; u++) {
			for (int k=0; k<my_req_num; k++) {
				obj += cs_len * ys[u][k] * xs[u][k];
			}
		}
	}

	model.add(IloMaximize(env, obj));
	obj.end();
}


// Constraint 9 (Priority): at most 1 request from all lower priority tasks 
// can precede a request sent from me
void BlockingAnalysis::add_prio_constraint_at_most_one_lp_request(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
																  map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
																  map<ResourceID, map<TaskID, CSData> > &x_data, 
																  IloRangeArray &cons, Task *task, TaskSet *taskset, IloEnv &env) {
	TaskID task_id = task->task_id_;
	vector<TaskID> &my_lp_tasks = taskset->lower_prio_tasks_[task_id];
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator it = x_vars.begin();
	for (; it != x_vars.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, IloNumVarArray> &interfere = it->second;
		IloExpr expr(env);
		for (map<TaskID, IloNumVarArray>::iterator iter=interfere.begin(); iter!=interfere.end(); iter++) {
			TaskID interfere_task_id = iter->first;
			if (find(my_lp_tasks.begin(), my_lp_tasks.end(), interfere_task_id) != my_lp_tasks.end()) {
				IloNumVarArray &lp_x = iter->second;
				unsigned int his_request_num = x_data[res_id][interfere_task_id].request_num;
				for (int i=0; i<his_request_num; i++) {
					expr += lp_x[i];
				}
			}
		}

		unsigned int my_request_num = task->my_resources_[res_id]->request_num;
		// considering the first processor (id=0)
		IloNumVarArray &my_y = my_y_vars[res_id][0];
		for (int i=0; i<my_request_num; i++) {
			expr -= my_y[i];
		}

		cons.add(expr <= 0);
		expr.end();
	}
}

// Calculate delay-per-request for priority-ordered locks with unordered tie break
double BlockingAnalysis::delay_per_request(Task *task, TaskSet *taskset, ResourceID rid, SpinlockType lock_type) {
	if (lock_type != PRIO_UNORDERED && lock_type != PRIO_FIFO) {
		cout << "Wrong lock type in calculation of delay-per-request !!!" << endl;
		return 0;
	}

	TaskID my_id = task->task_id_;
	map<TaskID, Task*> &tset = taskset->tasks_;
	vector<TaskID> &lp_tasks = taskset->lower_prio_tasks_[my_id];
	vector<TaskID> &hp_tasks = taskset->higher_prio_tasks_[my_id];
	vector<TaskID> &ep_tasks = taskset->equal_prio_tasks_[my_id];

	// List of interfering tasks w.r.t this resource
	vector<TaskID> *interfering_tasks = task->interferences_[rid];

	// Collect lower, higher, equal locking-priority tasks w.r.t this resource
	vector<TaskID> lp_tasks_for_resource;
	vector<TaskID> hp_tasks_for_resource;
	vector<TaskID> ep_tasks_for_resource;
	for (int i=0; i<interfering_tasks->size(); i++) {
		TaskID his_id = interfering_tasks->at(i);
		if (find(lp_tasks.begin(), lp_tasks.end(), his_id) != lp_tasks.end()) {
			lp_tasks_for_resource.push_back(his_id);
		} else if (find(hp_tasks.begin(), hp_tasks.end(), his_id) != hp_tasks.end()) {
			hp_tasks_for_resource.push_back(his_id);
		} else if (find(ep_tasks.begin(), ep_tasks.end(), his_id) != ep_tasks.end()) {
			ep_tasks_for_resource.push_back(his_id);
		} else {
			cout << "SOMETHING WRONG!!!" << endl;
		}
	}

	// Find the maximum length of a critical section among
	// lower locking-priority tasks (w.r.t this resource)
	double max_cslen_lp = 0;
	for (int i=0; i<lp_tasks_for_resource.size(); i++) {
		TaskID his_id = lp_tasks_for_resource[i];
		double his_cs_len = tset[his_id]->my_resources_[rid]->critical_section_len;
		if (his_cs_len > max_cslen_lp)
			max_cslen_lp = his_cs_len;
	}

	// Add 1 to avoid the case when max_cslen_lp equals 0
	double dpr = max_cslen_lp + 1;

	// Debugging variable to count the number of iterations
	int count = 0;

	bool converged = false;
	while (!converged) {
		count++;
		// Init tmp (interferences from LOWER locking-priority requests)
		double tmp = max_cslen_lp + 1;
		
		// Add interferences from HIGHER locking-priority requests
		for (int i=0; i<hp_tasks_for_resource.size(); i++) {
			TaskID his_id = hp_tasks_for_resource[i];
			double his_cs_len = tset[his_id]->my_resources_[rid]->critical_section_len;
			unsigned int his_req_num = tset[his_id]->my_resources_[rid]->request_num;
			tmp += njobs(tset[his_id], dpr) * his_req_num * his_cs_len;
		}

		// Add interferences from EQUAL locking-priority requests
		if (lock_type == PRIO_UNORDERED) {
			for (int i=0; i<ep_tasks_for_resource.size(); i++) {
				TaskID his_id = ep_tasks_for_resource[i];
				double his_cs_len = tset[his_id]->my_resources_[rid]->critical_section_len;
				unsigned int his_req_num = tset[his_id]->my_resources_[rid]->request_num;
				tmp += njobs(tset[his_id], dpr) * his_req_num * his_cs_len;
			}
			double my_cs_len = task->my_resources_[rid]->critical_section_len;
			unsigned int my_req_num = task->my_resources_[rid]->request_num;
			tmp += (my_req_num - 1)*my_cs_len;
		} else if (lock_type == PRIO_FIFO) {
			for (int i=0; i<ep_tasks_for_resource.size(); i++) {
				TaskID his_task_id = ep_tasks_for_resource[i];
				double his_cs_len = tset[his_task_id]->my_resources_[rid]->critical_section_len;
				unsigned int his_proc_num = tset[his_task_id]->proc_num_;
				unsigned int his_req_num = njobs(tset[his_task_id], task->response_time_) * tset[his_task_id]->my_resources_[rid]->request_num;
				tmp += min(his_proc_num, his_req_num) * his_cs_len;
			}
			unsigned int my_proc_num = task->proc_num_;
			unsigned int my_req_num = task->my_resources_[rid]->request_num;
			double my_cs_len = task->my_resources_[rid]->critical_section_len;
			tmp += min(my_proc_num-1, my_req_num-1) * my_cs_len;
		}

		// Set flag to true if dpr is converged
		if (tmp == dpr)
			converged = true;
		
		// Update value for dpr either way
		dpr = tmp;
	}

	return dpr;
}


// Constraint 10: bound contributions of higher and equal locking-priority tasks 
// for with unordered tie break
void BlockingAnalysis::add_prio_constraint_unordered_tiebreak(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
															  map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
															  map<ResourceID, map<TaskID, CSData> > &x_data, 
															  IloRangeArray &cons, Task *task, map<ResourceID, double> &dpr, TaskSet *taskset, IloEnv &env) {
	TaskID my_task_id = task->task_id_;
	vector<TaskID> &my_hp_tasks = taskset->higher_prio_tasks_[my_task_id];
	vector<TaskID> &my_ep_tasks = taskset->equal_prio_tasks_[my_task_id];
	map<TaskID, Task*> &tset = taskset->tasks_;

	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator it = x_vars.begin();
	for (; it != x_vars.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, IloNumVarArray> &interfere = it->second;
		for (map<TaskID, IloNumVarArray>::iterator iter = interfere.begin(); iter != interfere.end(); iter++) {
			TaskID his_task_id = iter->first;
			if (find(my_hp_tasks.begin(), my_hp_tasks.end(), his_task_id) != my_hp_tasks.end() ||
				find(my_ep_tasks.begin(), my_ep_tasks.end(), his_task_id) != my_ep_tasks.end()) {
				// Either higher or equal locking-priority tasks
				IloExpr lhs(env);
				unsigned int his_request_num = x_data[res_id][his_task_id].request_num;
				// Left-hand side of constraint 10
				for (int i=0; i<his_request_num; i++) {
					lhs += iter->second[i];
				}
				
				// Right-hand side of constraint 10
				IloExpr rhs(env);
				double my_dpr = dpr[res_id];
				int coeff = njobs(tset[his_task_id], my_dpr) * tset[his_task_id]->my_resources_[res_id]->request_num;
				unsigned int my_req_num = task->my_resources_[res_id]->request_num;
				IloNumVarArray &y_vars = my_y_vars[res_id][0];
				for (int i=0; i<my_req_num; i++) {
					rhs += coeff * y_vars[i];
				}
				
				// Add constraint
				cons.add(lhs - rhs <= 0);

				lhs.end();
				rhs.end();
			}
		}
	}
}

// Constraint for Priority-Unordered locks 
// This constraint dictates that the blocking caused by other processors
// of tau_i to B_{i,j} is 0 if there is no request sent from professor j(th)
void BlockingAnalysis::add_prio_constraint_unordered_tiebreak_self(map<ResourceID, vector<IloNumVarArray> > &my_y_vars,
																   map<ResourceID, vector<IloNumVarArray> > &my_x_vars,
																   IloRangeArray &cons, Task *task, IloEnv &env) {
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	unsigned int my_proc_num = task->proc_num_;
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID res_id = y_it->first;
		vector<IloNumVarArray> &y_vars = y_it->second;
		vector<IloNumVarArray> &x_vars = my_x_vars[res_id];
		int my_req_num = task->my_resources_[res_id]->request_num;

		IloExpr lhs(env);
		for (int i=1; i<my_proc_num; i++) {
			for (int j=0; j<my_req_num; j++) {
				lhs += y_vars[i][j]*x_vars[i][j];
			}
		}

		IloExpr rhs(env);
		for (int j=0; j<my_req_num; j++) {
			rhs += my_req_num * y_vars[0][j];
		}

		cons.add(lhs - rhs <= 0);
		lhs.end();
		rhs.end();
	}
}


// Constraint 11: contribution of higher locking-priority tasks for priority-ordered locks 
// with FIFO tie break
void BlockingAnalysis::add_prio_constraint_fifo_tiebreak_hp(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
															map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
															map<ResourceID, map<TaskID, CSData> > &x_data,
															IloRangeArray &cons, Task *task, map<ResourceID, double> &dpr, TaskSet *taskset, IloEnv &env) {
	TaskID my_task_id = task->task_id_;
	vector<TaskID> &my_hp_tasks = taskset->higher_prio_tasks_[my_task_id];
	map<TaskID, Task*> &tset = taskset->tasks_;
	
	map<ResourceID, map<TaskID, IloNumVarArray> >:: iterator it = x_vars.begin();
	for (; it != x_vars.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, IloNumVarArray> &interfere = it->second;
		for (map<TaskID, IloNumVarArray>::iterator iter = interfere.begin(); iter != interfere.end(); iter++) {
			TaskID his_task_id = iter->first;
			if (find(my_hp_tasks.begin(), my_hp_tasks.end(), his_task_id) != my_hp_tasks.end()) {
				// Interfering task is a higher locking-priority task
				IloExpr expr(env);
				unsigned int his_req_num = x_data[res_id][his_task_id].request_num;
				IloNumVarArray &his_x_vars = iter->second;
				for (int i=0; i<his_req_num; i++) {
					expr += his_x_vars[i];
				}

				// Right-hand side expression
				double res_dpr = dpr[res_id];
				int coeff = njobs(tset[his_task_id], res_dpr) * tset[his_task_id]->my_resources_[res_id]->request_num;
				unsigned int my_req_num = task->my_resources_[res_id]->request_num;
				IloNumVarArray &y_vars = my_y_vars[res_id][0];
				for (int i=0; i<my_req_num; i++) {
					expr -= coeff * y_vars[i];
				}

				// Add constraint
				cons.add(expr <= 0);
				expr.end();
			}
		}
	}
}


// Constraint 12: bound contribution of equal locking-priority tasks for 
// Priority-ordered locks with FIFO tie-break
void BlockingAnalysis::add_prio_constraint_fifo_tiebreak_ep(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars,
															map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
															map<ResourceID, map<TaskID, CSData> > &x_data,
															IloRangeArray &cons, Task *task, TaskSet *taskset, IloEnv &env) {
	vector<TaskID> &my_ep_tasks = taskset->equal_prio_tasks_[task->task_id_];
	map<TaskID, Task*> &tset = taskset->tasks_;
	
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator it = x_vars.begin();
	for (; it != x_vars.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, IloNumVarArray> &interfere = it->second;
		for (map<TaskID, IloNumVarArray>::iterator iter = interfere.begin(); iter != interfere.end(); iter++) {
			TaskID his_task_id = iter->first;
			if (find(my_ep_tasks.begin(), my_ep_tasks.end(), his_task_id) != my_ep_tasks.end()) {
				unsigned int his_req_num = x_data[res_id][his_task_id].request_num;
				IloExpr expr1(env);
				IloExpr expr2(env);
				// Left-hand side of contraint 12
				for (int i=0; i<his_req_num; i++) {
					expr1 += iter->second[i];
					expr2 += iter->second[i];
				}

				// Right-hand side of constraint 12
				expr1 -= his_req_num;
				cons.add(expr1 <= 0);
				expr1.end();
				unsigned int my_req_num = task->my_resources_[res_id]->request_num;
				int his_proc_num = tset[his_task_id]->proc_num_;
				for (int i=0; i<my_req_num; i++) {
					expr2 -= his_proc_num * my_y_vars[res_id][0][i];
				}
				cons.add(expr2 <= 0);
				expr2.end();
			}
		}
	}
}


// Constraint 13: bound contribution of requests from other processors of the same task
void BlockingAnalysis::add_prio_constraint_fifo_tiebreak_self(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
															  map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
															  IloRangeArray &cons, Task *task, IloEnv &env) {
	add_fifo_constraint_myself(my_y_vars, my_x_vars, cons, task, env);
}


// Calculate maximum blocking for FIFO locks using numerical approach
double BlockingAnalysis::max_blocking_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data) {
	unsigned int my_proc_num = task->proc_num_;
	double my_blocking = 0;
	map<ResourceID, map<TaskID, CSData> >::iterator it = x_data.begin();
	for (; it != x_data.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, CSData> &interfere = it->second;
		unsigned int my_req_num = task->my_resources_[res_id]->request_num;
		double cs_len = task->my_resources_[res_id]->critical_section_len;

		vector<double> blks;
		for (int i=0; i<=my_req_num; i++) {
			// Add up blocking inside the task
			double tmp = min((my_proc_num-1)*i, my_req_num-i)*cs_len;

			// Add up blocking from other tasks
			for (map<TaskID, CSData>::iterator inner_it=interfere.begin(); inner_it!=interfere.end(); inner_it++) {
				TaskID his_task_id = inner_it->first;
				unsigned int his_proc_num = inner_it->second.proc_num;
				unsigned int his_req_num = inner_it->second.request_num;
				double his_cs_len = inner_it->second.cs_len;
				tmp += min(his_proc_num*i, his_req_num)*his_cs_len;
			}

			blks.push_back(tmp);
		}
		my_blocking += *max_element(blks.begin(), blks.end());
	}
	return my_blocking;
}

// Sort the list of interference tasks by their critical section length 
// @input: interfere - map of all interference tasks and their cs data
// @output: results - sorted cs data of lower priority tasks by decreasing cs length
void BlockingAnalysis::sort_by_cslen(Task *task, TaskSet *taskset, map<TaskID, CSData> &interfere, vector<CSData> &results) {
	vector<TaskID> &my_lp_tasks = taskset->lower_prio_tasks_[task->task_id_];

	for (map<TaskID, CSData>::iterator it = interfere.begin(); it != interfere.end(); it++) {
		// Only consider lower priotiry tasks
		if (find(my_lp_tasks.begin(), my_lp_tasks.end(), it->first) == my_lp_tasks.end())
			continue;

		// Find the correct insert position
		double current_cslen = it->second.cs_len;
		vector<CSData>::iterator insert_pos = results.begin();
		for (; insert_pos!=results.end(); insert_pos++) {
			if (current_cslen < insert_pos->cs_len) {
				continue;
			} else {
				break;
			}
		}

		// Insert to the result list
		results.insert(insert_pos, it->second);
	}
}

// Maximum blocking for Priority-ordered locks using numerical method
double BlockingAnalysis::max_blocking_prio(Task *task, TaskSet *taskset, map<ResourceID, map<TaskID, CSData> > &x_data,
										   map<ResourceID, double> &dpr, SpinlockType lock_type) {
	map<TaskID, Task*> &tset = taskset->tasks_;
	vector<TaskID> &my_lp_tasks = taskset->lower_prio_tasks_[task->task_id_];
	vector<TaskID> &my_ep_tasks = taskset->equal_prio_tasks_[task->task_id_];
	vector<TaskID> &my_hp_tasks = taskset->higher_prio_tasks_[task->task_id_];
	unsigned int my_proc_num = task->proc_num_;
	double my_blocking = 0;

	map<ResourceID, map<TaskID, CSData> >::iterator it = x_data.begin();
	for (; it != x_data.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, CSData> &interfere = it->second;
		unsigned int my_req_num = task->my_resources_[res_id]->request_num;
		double cs_len = task->my_resources_[res_id]->critical_section_len;

		// Sort interference tasks' data by cs length (decreasing order)
		vector<CSData> sorted_by_cslen;
		sort_by_cslen(task, taskset, interfere, sorted_by_cslen);
		
		vector<double> blks;
		for (int i=0; i<=my_req_num; i++) {
			double tmp_blocking = 0;
			// Add up blocking inside the task
			if (lock_type == PRIO_UNORDERED) {
				if (i > 0)
					tmp_blocking = (my_req_num-i)*cs_len;
			} else if (lock_type == PRIO_FIFO) {
				tmp_blocking = min((my_proc_num-1)*i, my_req_num-i)*cs_len;
			}

			// Add up blocking from lower locking-priority tasks
			unsigned int remain_req = i;
			for (int i=0; i<sorted_by_cslen.size(); i++) {
				if (sorted_by_cslen[i].request_num >= remain_req) {
					tmp_blocking += remain_req * sorted_by_cslen[i].cs_len;
					break;
				} else {
					tmp_blocking += sorted_by_cslen[i].request_num * sorted_by_cslen[i].cs_len;
					remain_req -= sorted_by_cslen[i].request_num;
				}
			}

			if (lock_type == PRIO_UNORDERED) {
				map<TaskID, CSData>::iterator inter_it = interfere.begin();
				
				for (; inter_it != interfere.end(); inter_it++) {
					TaskID his_task_id = inter_it->first;

					if (find(my_ep_tasks.begin(), my_ep_tasks.end(), his_task_id) != my_ep_tasks.end() ||
						find(my_hp_tasks.begin(), my_hp_tasks.end(), his_task_id) != my_hp_tasks.end() ) {
						// Add blocking from equal & higher priority tasks
						unsigned int his_req_num = inter_it->second.request_num;
						double his_cs_len = inter_it->second.cs_len;
						double dpr_for_curr_res = dpr[res_id];
						unsigned int req_inside_dpr = njobs(tset[his_task_id], dpr_for_curr_res)*tset[his_task_id]->my_resources_[res_id]->request_num;
						tmp_blocking += min(req_inside_dpr*i, his_req_num) * his_cs_len;
					}
				}

			} else if (lock_type == PRIO_FIFO) {
				map<TaskID, CSData>::iterator inter_it = interfere.begin();

				for (; inter_it != interfere.end(); inter_it++) {
					TaskID his_task_id = inter_it->first;
					
					if (find(my_ep_tasks.begin(), my_ep_tasks.end(), his_task_id) != my_ep_tasks.end() ) {
						// Add blocking from equal priority tasks
						unsigned int his_req_num = inter_it->second.request_num;
						unsigned int his_proc_num = inter_it->second.proc_num;
						double his_cs_len = inter_it->second.cs_len;

						tmp_blocking += min(his_proc_num*i, his_req_num)*his_cs_len;
					} else if (find(my_hp_tasks.begin(), my_hp_tasks.end(), his_task_id) != my_hp_tasks.end() ) {
						// Add blocking from higher priority tasks
						unsigned int his_req_num = inter_it->second.request_num;
						double his_cs_len = inter_it->second.cs_len;
						double dpr_for_curr_res = dpr[res_id];
						unsigned int req_inside_dpr = njobs(tset[his_task_id], dpr_for_curr_res)*tset[his_task_id]->my_resources_[res_id]->request_num;
						tmp_blocking += min(req_inside_dpr*i, his_req_num) * his_cs_len;
					}
				}

			} else {
				cerr << "WRONG PRIORITY LOCK TYPE !!!" << endl;
				exit(-1);
			}

			blks.push_back(tmp_blocking);
		}

		my_blocking += *max_element(blks.begin(), blks.end());
	}

	return my_blocking;
}

// Using numerical method instead of optimization to bound blocking on a single processor of a task
void BlockingAnalysis::blocking_bound_single_proc(Task *task, TaskSet *taskset, unsigned int m, double *blocking, SpinlockType lock_type) {
	TaskID my_id = task->task_id_;
	TaskID r_i = task->response_time_;
	map<ResourceID, vector<TaskID>*> &interferences = task->interferences_;
	map<ResourceID, vector<TaskID>*>::iterator rit = interferences.begin();
	map<TaskID, Task*> &tset = taskset->tasks_;

	// Gather information from requests of tau_x which 
	// interfere with my requests
	map<ResourceID, map<TaskID, CSData> > x_data;
	for (; rit != interferences.end(); rit++) {
		ResourceID r_id = rit->first;
		vector<TaskID>* vec = rit->second;
		for (int i=0; i<vec->size(); i++) {
			TaskID t_id = vec->at(i);
			// Abort if this is me, but this is impossible
			if (t_id == my_id)
				continue;
			Task* tau_x = tset[t_id];
			unsigned int proc_num = tau_x->proc_num_;
			map<ResourceID, Resource*> &tau_x_resources = tau_x->my_resources_;
			unsigned int n_x_q = tau_x_resources[r_id]->request_num;
			double cs_len = tau_x_resources[r_id]->critical_section_len;
			unsigned int request_num = njobs(tau_x, r_i)*n_x_q;
			map<TaskID, CSData> &inner_map = x_data[r_id];
			CSData data = {proc_num, request_num, cs_len};
			inner_map.insert(std::pair<TaskID, CSData> (t_id, data));
		}
	}

	if (lock_type == FIFO) {
		// Maximum blocking for FIFO locks
		*blocking = max_blocking_fifo(task, x_data);
	} else if (lock_type == PRIO_UNORDERED || lock_type == PRIO_FIFO) {
		// Calculate delay-per-request if lock type is priority-ordered
		map<ResourceID, double> dpr;
		map<ResourceID, Resource*> &my_resources = task->my_resources_;
		for (map<ResourceID, Resource*>::iterator it = my_resources.begin(); it != my_resources.end(); it++) {
			ResourceID res_id = it->first;
			double delay = delay_per_request(task, taskset, res_id, lock_type);
			dpr[res_id] = delay;
		}
		
		// Maximum blocking for Priority locks
		*blocking = max_blocking_prio(task, taskset, x_data, dpr, lock_type);
	}
}


// Max blocking for a whole job with FIFO locks
double BlockingAnalysis::max_blocking_whole_job_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data) {
	TaskID my_task_id = task->task_id_;
	unsigned int my_proc_num = task->proc_num_;
	
	double blocking = 0;
	map<ResourceID, map<TaskID, CSData> >::iterator it = x_data.begin();
	for (; it != x_data.end(); it++) {
		ResourceID res_id = it->first;
		map<TaskID, CSData> &interfere = it->second;
		unsigned int my_req_num = task->my_resources_[res_id]->request_num;
		double my_cs_len = task->my_resources_[res_id]->critical_section_len;
		
		// Add up self blocking
		if (my_req_num < my_proc_num) {
			blocking += (my_req_num*(my_req_num-1)/2) * my_cs_len;
		} else {
			blocking += ((my_proc_num-1)*(my_req_num-my_proc_num) + my_proc_num*(my_proc_num-1)/2) * my_cs_len;
		}

		// Add up blocking from other tasks
		for (map<TaskID, CSData>::iterator inter_it=interfere.begin(); inter_it!=interfere.end(); inter_it++) {
			TaskID his_task_id = inter_it->first;
			unsigned int his_req_num = inter_it->second.request_num;
			double his_cs_len = inter_it->second.cs_len;
			unsigned int his_proc_num = inter_it->second.proc_num;

			blocking += min(my_proc_num*his_req_num, his_proc_num*my_req_num) * his_cs_len;
		}
	}

	return blocking;
}

// Blocking bound of a whole job
void BlockingAnalysis::blocking_bound_whole_job(Task *task, TaskSet *taskset, unsigned int m, double *blocking, SpinlockType lock_type) {
	TaskID my_id = task->task_id_;
	TaskID r_i = task->response_time_;
	map<ResourceID, vector<TaskID>*> &interferences = task->interferences_;
	map<ResourceID, vector<TaskID>*>::iterator rit = interferences.begin();
	map<TaskID, Task*> &tset = taskset->tasks_;

	// Gather information from requests of tau_x which 
	// interfere with my requests
	map<ResourceID, map<TaskID, CSData> > x_data;
	for (; rit != interferences.end(); rit++) {
		ResourceID r_id = rit->first;
		vector<TaskID>* vec = rit->second;
		for (int i=0; i<vec->size(); i++) {
			TaskID t_id = vec->at(i);
			// Abort if this is me, but this is impossible
			if (t_id == my_id)
				continue;
			Task* tau_x = tset[t_id];
			unsigned int proc_num = tau_x->proc_num_;
			map<ResourceID, Resource*> &tau_x_resources = tau_x->my_resources_;
			unsigned int n_x_q = tau_x_resources[r_id]->request_num;
			double cs_len = tau_x_resources[r_id]->critical_section_len;
			unsigned int request_num = njobs(tau_x, r_i)*n_x_q;
			map<TaskID, CSData> &inner_map = x_data[r_id];
			CSData data = {proc_num, request_num, cs_len};
			inner_map.insert(std::pair<TaskID, CSData> (t_id, data));
		}
	}

	if (lock_type == FIFO) {
		*blocking = max_blocking_whole_job_fifo(task, x_data);
	} else if (lock_type == PRIO_UNORDERED || lock_type == PRIO_FIFO) {
		// Not support yet
	}
}

// Schedulability test using the following procedure:
// In each iteration, update blocking bound, then update number of processors 
// allocated for tasks. For each task, if the new number of processors allocated 
// to it decreases, keep it as the previous iteration. Otherwise, if it increases,
// update it. The test terminates when the total number of processors > m (unschedulable)
// OR processor allocations for all tasks do not change (schedulable)

// Number of times a taskset is unschedulable caused by a task in the taskset
// has inflated critical path length larger than D (L+Bi > D)
//int large_critical_path_length_fail_num=0;
//int different_num_opt_numerical = 0;
//int numerical_larger_opt_num = 0;
//int opt_larger_numerical_num = 0;
//int equal_num_opt_numerical = 0;
bool BlockingAnalysis::is_schedulable(TaskSet *tset, unsigned int m, SpinlockType lock_type) {
	map<TaskID, Task*> &taskset = tset->tasks_;
	map<TaskID, Task*>::iterator it = taskset.begin();

	// For all tasks, set their response times to relative deadline
	for (; it != taskset.end(); it++) {
		it->second->response_time_ = it->second->period_;
	}

	// Copy values for new_proc_num_ fields
	for (it = taskset.begin(); it != taskset.end(); it++) {
		it->second->new_proc_num_ = it->second->proc_num_;
	}

	//	int iter_num = 0;
	while (true) {
		//		iter_num++;
		unsigned int total_proc_num = 0;
		//		SCIP_Real blocking;
		double blocking;
		
		// Copy values of newProcNum back to procNum for each task
		for (it = taskset.begin(); it != taskset.end(); it++) {
			if (it->second->converged_ == false)
				it->second->proc_num_ = it->second->new_proc_num_;
		}
		
		// Calculate new processors allocations for tasks
		for (it = taskset.begin(); it != taskset.end(); it++) {
			Task * task = it->second;
			//			task_analysis(task, tset, m, &blocking, lock_type);

			// Blocking bound on a single processor using numerical method
			blocking_bound_single_proc(task, tset, m, &blocking, lock_type);
			//			if ( fabs(blocking - blk_numerical) > 0.0001) {
			//				different_num_opt_numerical++;
			//				cout << "Blocking by SCIP: " << blocking << "; Blocking by numerical: " << blk_numerical << endl;
			//				if (fabs(blocking - blk_numerical) >= 0.5 && blocking > blk_numerical) {
			//					opt_larger_numerical_num++;
			//				} else if (fabs(blocking - blk_numerical) >= 0.5 && blocking < blk_numerical) {
			//					numerical_larger_opt_num++;
			//				}
			//			} else {
			//				equal_num_opt_numerical++;
			//			}
			
			/* If D < L + B1, this task cannot be schedulable */
			if (task->period_ < task->span_ + blocking) {
				//				large_critical_path_length_fail_num++;
				return false;
			}

			/*
			// Blocking bound of the whole job
			double blocking_whole_job;
			if (lock_type == FIFO) {
				blocking_bound_whole_job(task, tset, m, &blocking_whole_job, lock_type);

				cout << "Whole blocking: " << blocking_whole_job << "; Old whole blocking: " << task->procNum*blocking << endl;
				if (blocking_whole_job > task->procNum*blocking) {
					cout << "New whole job's blocking is worse!" << endl;
				}
			}

			// New number of processors needed
			unsigned int proc_num;
			if (lock_type == FIFO) {
				proc_num = alloc_proc(task->C, task->L, task->T, blocking_whole_job, blocking);
			} else {
				proc_num = alloc_proc(task->C, task->L, task->T, task->procNum*blocking, blocking);
			}
			*/
			unsigned int proc_num = alloc_proc(task->wcet_, task->span_, task->period_, task->proc_num_*blocking, blocking);
			
			
			// If the new processors allocated is larger than the old one, update it,
			// otherwise, do not update it
			if (proc_num > task->proc_num_) {
				task->new_proc_num_ = proc_num;
				task->converged_ = false;
			} else {
				task->converged_ = true;
			}
			
			// Update task's blocking bound, total number of processors allocated so far
			task->blocking_on_single_ = blocking;
			total_proc_num += task->new_proc_num_;
		}

		//		cout << "Iter " << iter_num << "; total #processors: " << total_proc_num << endl;
		
		// If total processors used is larger than m, the taskset is unschedulable
		if (total_proc_num > m) {
			//			cout << "Total processors needed > m" << endl;
			return false;
		}

		// If processor allocations for all tasks do not change from the previous iteration
		// then the taskset is considered as converged
		bool global_converged = true;
		for (it = taskset.begin(); it != taskset.end(); it++) {
			global_converged &= it->second->converged_;
		}

		if (global_converged == true) {
			//			cout << "Global converged reaches" << endl;
			return true;
		}
	}
	
}


// Return number of tasks tau_x in interval length t
unsigned int BlockingAnalysis::njobs(Task* tau_x, double t) {
	double tau_x_period = tau_x->period_;
	double tau_x_response_time = tau_x->response_time_;	
	return ceil( (t + tau_x_response_time) / tau_x_period);
}

#if 0
// Try if this will be less pessimistic
unsigned int BlockingAnalysis::njobs(Task *tau_x, double t) {
	return (floor(t/tau_x->period_) + 1);
}
#endif
