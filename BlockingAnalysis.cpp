#include <math.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "BlockingAnalysis.h"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
	
ILOSTLBEGIN

/* Update number of processors allocated */
unsigned int alloc_proc(double C, double L, double D, double Bn, double B1) {
	if (D < L+B1)
		cout << "Denominator of equation for n_i is negative !!!" << endl;
	else if (D == L+B1)
		cout << "Denominator of equation for n_i is zero !!!" << endl;

	cout << "New processor allocation: " << ceil(((C+Bn) - (L+B1))/(D - (L+B1))) << endl;
	return ceil(((C+Bn) - (L+B1))/(D - (L+B1)));
}

/* Update bound of response time using blocking time, #processors */
double response_time(double C, double L, unsigned int n, double Bn, double B1) {
	double respTime = (C+Bn-L-B1)/n + L + B1;
	return ceil(respTime); // in microsecond
}

/* Initialize blocking time, #processors, response time 
 * Also, update "interferences" for each task in task set
 * Return: true, if taskset are OK after initiating
 *         false, otherwise
 * (OK means total processors allocated to tasks < m AND 
 *  there is no task with Response time > Deadline)
 */
bool init_iteration(TaskSet* taskset, unsigned int m) {
	map<TaskID, Task*> &tset = taskset->tasks;
	map<TaskID, Task*>::iterator it = tset.begin();
	for (; it != tset.end(); it++) {
		Task* task = it->second;
		task->B1 = 0;
		task->procNum = alloc_proc(task->C, task->L, task->T, 0, 0);
		task->R = response_time(task->C, task->L, task->procNum, 0, 0);
		task->converged = false;
		
		// Update interferences for this task
		map<ResourceID, vector<TaskID>*> &interferences = task->interferences;
		map<ResourceID, Resource*> &myResources = task->myResources;
		map<TaskID, Task*>::iterator taskIter;
		map<ResourceID, Resource*>::iterator resIter;

		for (resIter=myResources.begin(); resIter!=myResources.end(); resIter++) {
			ResourceID rid = resIter->first;
			interferences.insert(std::pair<ResourceID, vector<TaskID>*> 
								 (rid, new vector<TaskID>()));
			for (taskIter=tset.begin(); taskIter!=tset.end(); taskIter++) {
				// Abort if this is me
				if (taskIter->second->taskID == task->taskID)
					continue;

				// Add tasks access to the same resource
				if (taskIter->second->myResources.find(rid) != 
					taskIter->second->myResources.end()) {
					interferences.find(rid)->second->push_back(taskIter->first);
				}
			}
		}
	}


#define _INIT_DEBUG_
#define _INIT_VERBOSE_

#ifdef _INIT_DEBUG_
	/* Debug: dump tasks information */
	unsigned int total_proc = 0;
	bool isOk = true;
int number_tasks_fail = 0;
	for (it=tset.begin(); it!=tset.end(); it++) {
#ifdef _INIT_VERBOSE_
		cout << "Task ID: " << it->second->taskID << endl;
		cout << "Parameters (T,L,U,C): (" << it->second->T << "," <<
			it->second->L << "," << it->second->U << "," << 
			it->second->C << ")" << endl;
		cout << "#processors: " << it->second->procNum << endl;
		cout << "Response time: " << it->second->R << endl;
		cout << "Response time > Deadline: " <<
			((it->second->R > it->second->T)? "true" : "false") << endl;
#endif // _INIT_VERBOSE_
		
		total_proc += it->second->procNum;
		if (it->second->R > it->second->T) {
			number_tasks_fail++;
			isOk = false;
		}
		map<ResourceID, Resource*> &res = it->second->myResources;
		map<ResourceID, Resource*>::iterator resIt = res.begin();
#ifdef _INIT_VERBOSE_
		cout << "Resource ID list: (";
		for (; resIt!=res.end(); resIt++) {
			cout << resIt->first << " ";
		}
		cout << ")" << endl;

		map<ResourceID, vector<TaskID>*> &interfe = it->second->interferences;
		map<ResourceID, vector<TaskID>*>::iterator intIt = interfe.begin();
		cout << "Interference List ============= " << endl;
		for (; intIt!=interfe.end(); intIt++) {
			cout << "   For resource ID " << intIt->first << ": ";
			vector<TaskID>* vec = intIt->second;
			for (int i=0; i<vec->size(); i++) {
				cout << vec->at(i) << " ";
			}
			cout << endl;
		}
		cout << endl;
#endif // _INIT_VERBOSE_
	}
	cout << "Total #processors: " << total_proc << endl;
	if (total_proc > m || isOk == false) {
		cout << "#tasks fail: " << number_tasks_fail << "; Taskset fail !!!" << endl;
		return false;
	}
	return true;
#endif // _INIT_DEBUG_

	/* In case of no debugging, it just returns true 
	 * But it does not necessarily mean the task set is OK
	 */
	return true;
}

/* Structure to store information of requests 
 * from tau_x to l_q which intefere with tau_i
 * - requestNum: #requests of tau_x to l_q while tau_i pending
 * - csLen: CS length of request from tau_x to l_q (microsecond)
 */
typedef struct {
	unsigned int procNum;
	unsigned int requestNum;
	double csLen;
} CSData;

void add_generic_contraint_sum_of_y_for_each_request(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
													 IloRangeArray &cons, Task *task, IloEnv &env);
void add_generic_constraint_rule_out_mismatched_x_and_y(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
														map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
														IloRangeArray &cons, Task *task, IloEnv &env);
void add_fifo_constraint_other_tasks(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
									 map<ResourceID, map<TaskID, CSData> > &x_data, 
									 map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
									 IloRangeArray &cons, Task *task, IloEnv &env);
void add_fifo_constraint_myself(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
								map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
								IloRangeArray &cons, Task *task, IloEnv &env);
void add_objective_function(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
							map<ResourceID, map<TaskID, CSData> > &x_data, 
							map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
							map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
							Task *task, IloModel &model, IloEnv &env);

void call_optimizer_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data, SCIP_Real *max_blocking);

/* Update task's blocking time, #processors, and response time 
 * based on results from the previous iteration.
 */
void task_analysis(Task* task, TaskSet* taskset, unsigned int m, SCIP_Real *max_blocking) {
	TaskID myId = task->taskID;
	TaskID r_i = task->R;
	map<ResourceID, vector<TaskID>*> &interferences = task->interferences;
	map<ResourceID, vector<TaskID>*>::iterator rit = interferences.begin();
	map<TaskID, Task*> &tset = taskset->tasks;

	/* Gather information from requests of tau_x which 
	 * interfere with my requests
	 */
	map<ResourceID, map<TaskID, CSData> > x_data;
	for (; rit != interferences.end(); rit++) {
		ResourceID rId = rit->first;
		vector<TaskID>* vec = rit->second;
		for (int i=0; i<vec->size(); i++) {
			TaskID tId = vec->at(i);
			/* Abort if this is me, but this is impossible */
			if (tId == myId)
				continue;
			Task* tau_x = tset[tId];
			unsigned int procNum = tau_x->procNum;
			map<ResourceID, Resource*> &tau_x_resources = tau_x->myResources;
			unsigned int N_x_q = tau_x_resources[rId]->requestNum;
			double csLen = tau_x_resources[rId]->CSLength;
			unsigned int request_num = njobs(tau_x, r_i)*N_x_q;
			map<TaskID, CSData> &inner_map = x_data[rId];
			CSData data = {procNum, request_num, csLen};
			inner_map.insert(std::pair<TaskID, CSData> (tId, data));
		}
	}

	/* Print debug information */
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

	// Call the optimizer to solve the maximum blocking time
	SCIP_Real max_blk;
	call_optimizer_fifo(task, x_data, &max_blk);
	*max_blocking = max_blk;
}

SCIP_RETCODE optimize(string fname, SCIP_Real *obj_value);

void call_optimizer_fifo(Task *task, map<ResourceID, map<TaskID, CSData> > &x_data, SCIP_Real *blocking) {
	/* Solve the maximization of blocking */
	IloEnv env;
	try {
		//		cout << "Start modeling..." << endl;
		IloModel model(env);
		IloRangeArray cons(env);
		/* Variables table of requests from other tasks */
		map<ResourceID, map<TaskID, IloNumVarArray> > x_vars;
		map<ResourceID, map<TaskID, CSData> >::iterator iter = x_data.begin();

		//		cout << "TaskID: " << task->taskID << ", #Resources: " << x_data.size() << endl;

		//		cout << "Start create variables for interfering tasks" << endl;
		for (; iter != x_data.end(); iter++) {
			ResourceID resId = iter->first;
			//			cout << "Resource ID: " << resId << endl;
			/* Get an inner map for this resource ID */
			map<TaskID, IloNumVarArray> &innerMap = x_vars[resId];
			map<TaskID, CSData> &taskData = iter->second;
			map<TaskID, CSData>::iterator innerIter = taskData.begin();
			for (; innerIter != taskData.end(); innerIter++) {
				/* Each other task has a corresponding IloNumVarArray */
				TaskID taskId = innerIter->first;
				/* Get number of interfering requests from this task 
				 * to this resource
				 */
				unsigned int varNum = innerIter->second.requestNum;
				IloNumVarArray x(env);
				/* All x variables are in range [0,1] */
				for (int i=0; i<varNum; i++) {
					// Trick: Use a very small number (instead of 0) for lower bound of x vars
					//					x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
					// Use 0 as lower bound of x vars
					x.add(IloNumVar(env, 0, 1.0, ILOFLOAT));
				}
				innerMap.insert(std::pair<TaskID, IloNumVarArray>(taskId, x));
			} /* Task */
		} /* Resource */
		//		cout << "Stop creating variables for interfering tasks" << endl;
		
		/* Populate data for my (tau_i's) y & x variables 
		 * For my_y_vars: each element of a vector 
		 * corresponds to a row of matrix Y (n_i * N_i_q)
		 * For my_x_vars: each element of a vector 
		 * corresponds to a column of matrix X (N_i_q * n_i)
		 */
		unsigned int my_proc_num = task->procNum;
		map<ResourceID, vector<IloNumVarArray> > my_y_vars;
		map<ResourceID, vector<IloNumVarArray> > my_x_vars;
		map<ResourceID, Resource*> &myRes = task->myResources;
		map<ResourceID, Resource*>::iterator my_rit = myRes.begin();
		
		//		cout << "Start populating my variables" << endl;
		for (; my_rit != myRes.end(); my_rit++) {
			ResourceID rId = my_rit->first;
			//			cout << "Resource ID: " << rId << endl;
			unsigned int my_req_num = my_rit->second->requestNum;
			vector<IloNumVarArray> &ys = my_y_vars[rId];
			vector<IloNumVarArray> &xs = my_x_vars[rId];
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
		//		cout << "Stop populating my variables" << endl;

		/* Constraint 3: sum of each y column is <= 1 */
		add_generic_contraint_sum_of_y_for_each_request(my_y_vars, cons, task, env);

		/* Constraint 5: mismatched product of y&x is zero */
		add_generic_constraint_rule_out_mismatched_x_and_y(my_y_vars, my_x_vars, cons, task, env);

		/* Constraint 8: at most 1 request from each other processor of another task block me */
		add_fifo_constraint_other_tasks(x_vars, x_data, my_y_vars, cons, task, env);

		/* Constraint 9: doing the same for other processors in the same task */
		add_fifo_constraint_myself(my_y_vars, my_x_vars, cons, task, env);
		
		/* Add constraints to Cplex model */
		model.add(cons);
		
		/* Formulate the objective function & add obj function to Cplex model */
		add_objective_function(x_vars, x_data, my_y_vars, my_x_vars, task, model, env);

		IloCplex cplex(model);

		// disable output to stdout
		cplex.setOut(env.getNullStream());

		// create file name
		TaskID task_id = task->taskID;
		string file_name;
		stringstream ss;
		ss << "tmp/fifo_" << task_id << ".lp";
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

		/*
		if (!cplex.solve()) {
			env.error() << "Failed to optimize" << endl;
			throw(-1);
		} else {
			cout << "SOLVE SUCCESSFULLY" << endl;
		}

		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;

		/* Print out variables' values * /
		for (map<ResourceID, map<TaskID, IloNumVarArray> >::iterator it = x_vars.begin(); 
			 it != x_vars.end(); it++) {
			ResourceID rid = it->first;
			map<TaskID, IloNumVarArray> &inner_map = it->second;
			env.out() << "Inteferences from other tasks for resource ID " << rid << ": " << endl;

			for (map<TaskID, IloNumVarArray>::iterator inner_it = inner_map.begin();
				 inner_it != inner_map.end(); inner_it++) {
				TaskID tid = inner_it->first;
				IloNumVarArray &var_array = inner_it->second;
				IloNumArray val_array(env);
				cplex.getValues(val_array, var_array);
				env.out() << "Task ID " << tid << ": " << val_array << endl;
			}
		}

		for (map<ResourceID, vector<IloNumVarArray> >::iterator it = my_y_vars.begin();
			 it != my_y_vars.end(); it++) {
			ResourceID rid = it->first;
			vector<IloNumVarArray> &y_matrix = it->second;
			vector<IloNumVarArray> &x_matrix = my_x_vars[rid];
			env.out() << "Variables for resource ID: " << rid << endl;

			for (int i=0; i<y_matrix.size(); i++) {
				IloNumVarArray &y_row = y_matrix[i];
				IloNumVarArray &x_col = x_matrix[i];
				IloNumArray y_vals(env);
				IloNumArray x_vals(env);
				cplex.getValues(y_vals, y_row);
				cplex.getValues(x_vals, x_col);
				env.out() << "Processor #" << i << ": Ys = " << y_vals << "; Xs = " << x_vals << endl;
			}
		}
		*/

		/* Export the model to this file */
		//		cplex.exportModel("fifo.lp");

	} catch (IloException &e) {
		cerr << "Ilog concert exception: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception: " << endl;
	}

	cout << endl;
	env.end();
}

SCIP_RETCODE optimize(string fname, SCIP_Real *obj_value) {
	const char* file_name= fname.c_str();
	SCIP *scip = NULL;
	SCIP_CALL( SCIPcreate(&scip) );
	
	// load default plugins
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
	
	// disable output to stdout
	SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );
	
	SCIP_CALL( SCIPreadProb(scip, file_name, NULL) );
	
	//	std::cout << "Solving problem" << std::endl;
	SCIP_CALL( SCIPsolve(scip) );
	
	// get best solution
	SCIP_Real obj_val = SCIPgetPrimalbound(scip);
	//	std::cout << "Best solution for the problem: " << std::endl;
	std::cout << "Objective value: " << obj_val << std::endl;

	*obj_value = obj_val;
	
	SCIP_CALL( SCIPfree(&scip) );
	
	BMScheckEmptyMemory();
	return SCIP_OKAY;
}

/* Constraint 3: for each request to a resource, sum 
 * of corresponding y variables is at most 1
 * Constraint 2 then follows
 */
void add_generic_contraint_sum_of_y_for_each_request(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
													 IloRangeArray &cons, Task *task, IloEnv &env) {	
	//		cout << "Start imposing constraint 3" << endl;
	unsigned int my_proc_num = task->procNum;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	map<ResourceID, Resource*> &myRes = task->myResources;
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID rId = y_it->first;			
		vector<IloNumVarArray> &var_arrays = y_it->second;
		unsigned int req_num = myRes[rId]->requestNum;
		for (int k=0; k<req_num; k++) {
			IloExpr expr(env);
			for (int u=0; u<my_proc_num; u++) {
				expr += var_arrays[u][k];
			}
			//			cons.add(1 <= expr <= 1);
			cons.add(expr <= 1);
			expr.end(); /* Important: must free it */
		}
	}
	//		cout << "End imposing constraint 3" << endl;
}

/* Constraint 5: products of mismatched y and x variables 
 * are all zero
 */
void add_generic_constraint_rule_out_mismatched_x_and_y(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
														map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
														IloRangeArray &cons, Task *task, IloEnv &env) {
	//		cout << "Start imposing constraints 5 and 6" << endl;
	unsigned int my_proc_num = task->procNum;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	map<ResourceID, Resource*> &myRes = task->myResources;
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID rId = y_it->first;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[rId];
		unsigned int req_num = myRes[rId]->requestNum;
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
		
		/* Constraint 6: requests from me cannot block me 
		 * (processor under consideration is the first 
		 * processor of tau_i)
		 */
		IloExpr expr(env);
		for (int k=0; k<req_num; k++) {
			expr += ys[0][k] * xs[0][k];
		}
		cons.add(expr <= 0);
		expr.end();
	}
	//		cout << "Stop imposing constraints 5 and 6" << endl;
}


/* Constraint 8 (FIFO): at most 1 request from each other
 * processor of another task can block a single request from me
 */
void add_fifo_constraint_other_tasks(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
									 map<ResourceID, map<TaskID, CSData> > &x_data, 
									 map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
									 IloRangeArray &cons, Task *task, IloEnv &env) {
	//		cout << "Start imposing constraint 8 (for other tasks)" << endl;
	map<ResourceID, Resource*> &myRes = task->myResources;
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator inter_it = x_vars.begin();
	for (; inter_it != x_vars.end(); inter_it++) {
		ResourceID rId = inter_it->first;
		map<TaskID, IloNumVarArray> &others = inter_it->second;
		map<TaskID, IloNumVarArray>::iterator others_it = others.begin();
		for (; others_it != others.end(); others_it++) {
			TaskID tId = others_it->first;
			IloNumVarArray &x = others_it->second;
			IloExpr expr(env);
			unsigned int his_request_num = x_data[rId][tId].requestNum;
			int his_proc_num = x_data[rId][tId].procNum;
			for (int k=0; k<his_request_num; k++) {
				expr += x[k];
			}
			vector<IloNumVarArray> &my_ys = my_y_vars[rId];
			unsigned int my_request_num = myRes[rId]->requestNum;
			for (int i=0; i<my_request_num; i++) {
				expr -= his_proc_num * my_ys[0][i];
			}
			cons.add(expr <= 0);
			expr.end();
		}
	}
	//	cout  << "Stop imposing constraint 8 (for other tasks)" << endl;
}

/*
 * Constraint 8 (FIFO): for interference within myself
 */
void add_fifo_constraint_myself(map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
								map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
								IloRangeArray &cons, Task *task, IloEnv &env) {
	//	cout << "Start imposing constraint 8 (myself)" << endl;
	unsigned int my_proc_num = task->procNum;
	map<ResourceID, Resource*> &myRes = task->myResources;	
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID rId = y_it->first;
		unsigned int my_request_num = myRes[rId]->requestNum;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[rId];
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
	//		cout << "Stop imposing constraint 8 (myself)" << endl;
}

/*
 * Build objective function for the blocking time
 */
void add_objective_function(map<ResourceID, map<TaskID, IloNumVarArray> > &x_vars, 
							map<ResourceID, map<TaskID, CSData> > &x_data, 
							map<ResourceID, vector<IloNumVarArray> > &my_y_vars, 
							map<ResourceID, vector<IloNumVarArray> > &my_x_vars, 
							Task *task, IloModel &model, IloEnv &env) {
	//	cout << "Start building objective function" << endl;
	IloExpr obj(env);
	map<ResourceID, map<TaskID, IloNumVarArray> >::iterator inter_it = x_vars.begin();
	for (; inter_it != x_vars.end(); inter_it++) {
		ResourceID rId = inter_it->first;
		map<TaskID, IloNumVarArray> &others = inter_it->second;
		map<TaskID, IloNumVarArray>::iterator others_it = others.begin();
		for (; others_it != others.end(); others_it++) {
			TaskID tId = others_it->first;
			double csLen = x_data[rId][tId].csLen;
			unsigned int varNum = x_data[rId][tId].requestNum;
			for (int i=0; i<varNum; i++) {
				obj += csLen * others[tId][i];
			}
		}
	}
	
	unsigned int my_proc_num = task->procNum;
	map<ResourceID, Resource*> &myRes = task->myResources;
	map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
	for (; y_it != my_y_vars.end(); y_it++) {
		ResourceID rId = y_it->first;
		vector<IloNumVarArray> &ys = y_it->second;
		vector<IloNumVarArray> &xs = my_x_vars[rId];
		unsigned int my_req_num = myRes[rId]->requestNum;
		double csLen = myRes[rId]->CSLength;
		
		for (int u=0; u<my_proc_num; u++) {
			for (int k=0; k<my_req_num; k++) {
				obj += csLen * ys[u][k] * xs[u][k];
			}
		}
	}

	model.add(IloMaximize(env, obj));
	obj.end();
	//		cout << "Finish building objective function" << endl;
}


#define _SENSITIVITY_ 100
/**
 * Doing schelability analysis for a task set
 * Return true if it is schedulable, false otherwise
 */
bool blocking_analysis(TaskSet* tset, unsigned int m) {
	map<TaskID, Task*> &taskset = tset->tasks;
	map<TaskID, Task*>::iterator it = taskset.begin();
	
	int iter_num = 0;
	SCIP_Real blocking;
	while (true) {
		// Reset task set iterator
		it = taskset.begin();
		cout << "Iteration #" << iter_num << "XXXXXXXXXXXXXXXXXX" << endl;
		for (; it != taskset.end(); it++) {			
			Task* task = it->second;

			// Debug: print task's parameters
			cout << "Task " << task->taskID << ": Response time " << task->R 
				 << "; convergence: " << task->converged << endl;
			// End debugging

			if (task->converged == true)
				continue;
			
			task_analysis(task, tset, m, &blocking);
			double total_blk = blocking * task->procNum;
			unsigned int processor_num = alloc_proc(task->C, task->L, task->T, total_blk, blocking);
			double resp_time = response_time(task->C, task->L, processor_num, total_blk, blocking);
			
			// Debug: print total blocking, new processor number, new response time
			cout << "New blocking: " << blocking << "; new #proc: " << processor_num 
				 << "; new response time: " << resp_time << endl;
			// End debugging

			// if the new response time is larger than deadline, task set is unschedulable
			if (resp_time > task->T) {
				// Debug:
				cout << "New response time " << resp_time << " > period " << task->T 
					 << " ==> unschedulable" << endl;
				// End debugging
				return false;
			}

			double old_response_time = task->R;
			
			// if the response time does not change, the task is converged
			if ( fabs(resp_time - old_response_time) < _SENSITIVITY_ ) 
				task->converged = true;

			// update task information
			task->B1 = blocking;
			task->procNum = processor_num;
			task->R = resp_time;
		}

		unsigned int total_processors_allocated = 0;
		it = taskset.begin();
		for (; it != taskset.end(); it++) {
			total_processors_allocated += it->second->procNum;
		}
		
		// Debug:
		cout << "Total number of processors allocated: " << total_processors_allocated << endl;
		// End debugging
		
		// if the total number of processors allocated so far is more than m, task set is unschedulable
		if (total_processors_allocated > m) {
			// Debug:
			cout << "Not enough processors: m = " << m << " < " << total_processors_allocated << endl;
			// End debugging
			return false;
		}

		bool all_converged = true;
		it = taskset.begin();
		for (; it != taskset.end(); it++) {
			all_converged &= it->second->converged;
			if ( !all_converged )
				break;
		}
		
		if (all_converged)
			return true;

		// Debug: to keep track number of iterations
		iter_num++;
	}

}


/* For each task, fix its response time to its relative 
 * deadline and calculate the number of processors allocated.
 * If the total number of processors allocated to the taskset 
 * is less than m, than it is schedulable (only sufficient condition).
 */
bool is_schedulable(TaskSet *tset, unsigned int m) {
	map<TaskID, Task*> &taskset = tset->tasks;
	map<TaskID, Task*>::iterator it = taskset.begin();

	/* For all tasks, set their response times to relative deadline*/
	for (; it != taskset.end(); it++) {
		it->second->R = it->second->T;
	}

	unsigned int total_proc_num = 0;
	SCIP_Real blocking;
	it = taskset.begin();
	for (; it != taskset.end(); it++) {
		Task * task = it->second;
		task_analysis(task, tset, m, &blocking);
		
		/* If D <= L + B1, this task cannot be schedulable */
		if (task->T <= task->L + blocking) {
			return false;
		}
		unsigned int proc_num = alloc_proc(task->C, task->L, task->T, 
										   task->procNum*blocking, blocking);
		
		/* Update task's blocking bound, number of processors */
		task->B1 = blocking;
		task->procNum = proc_num;
		total_proc_num += proc_num;
	}
	
	cout << "New total processors allocated: " << total_proc_num << endl;
	
	if (total_proc_num <= m)
		return true;
	return false;
}


/* Helper funtions */

/* Return number of tasks tau_x in interval length t */
unsigned int njobs(Task* tau_x, double t) {
	double tau_x_period = tau_x->T;
	double tau_x_responsetime = tau_x->R;
	
	unsigned int ret = ceil((t+tau_x_responsetime)/tau_x_period);
	//	cout << "njobs() returns: " << ret << endl;
	return ret;
}
