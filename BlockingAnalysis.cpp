#include <math.h>
#include <iostream>
#include <ilcplex/ilocplex.h>
#include "BlockingAnalysis.h"

ILOSTLBEGIN

/* Update number of processors allocated */
unsigned int alloc_proc(double C, double L, double D, double Bn, double B1) {
	return ceil(((C+Bn) - (L+B1))/(D - (L+B1)));
}

/* Update bound of response time using blocking time, #processors */
double response_time(double C, double L, unsigned int n, double Bn, double B1) {
	double respTime = (C+Bn)/n + L + B1;
	return ceil(respTime); // in microsecond
}

/* Initialize blocking time, #processors, response time 
 * Also, update "interferences" for each task in task set
 * Return: true, if task set are OK after initiating
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


	//#define _INIT_DEBUG_
	//#define _INIT_VERBOSE_
#ifdef _INIT_DEBUG_
	/* Debug: dump tasks information */
	unsigned int total_proc = 0;
	bool isOk = true;
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
		if (it->second->R > it->second->T)
			isOk = false;
		map<ResourceID, Resource*> &res = it->second->myResources;
		map<ResourceID, Resource*>::iterator resIt = res.begin();
#ifdef _INIT_VERBOSE_
		cout << "Resource ID list: (";
		for (; resIt!=res.end(); resIt++) {
			cout << resIt->first << " ";
		}
		cout << ")" << endl;
#endif // _INIT_VERBOSE_

		map<ResourceID, vector<TaskID>*> &interfe = it->second->interferences;
		map<ResourceID, vector<TaskID>*>::iterator intIt = interfe.begin();
#ifdef _INIT_VERBOSE_
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
		cout << "Task fail !!!" << endl;
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

/* Update task's blocking time, #processors, and response time 
 * based on results from the previous iteration.
 */
void task_analysis(Task* task, TaskSet* taskset, unsigned int m) {
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

	/* Solve the maximization of blocking */
	IloEnv env;
	try {
		//		cout << "Start modeling..." << endl;
		IloModel model(env);
		IloRangeArray cons(env);
		/* Variables table of requests from other tasks */
		map<ResourceID, map<TaskID, IloNumVarArray> > x_vars;
		map<ResourceID, map<TaskID, CSData> >::iterator iter = x_data.begin();

		//		cout << "Start create variables for interfering tasks" << endl;
		for (; iter != x_data.end(); iter++) {
			ResourceID resId = iter->first;
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
					x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
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
			unsigned int my_req_num = my_rit->second->requestNum;
			vector<IloNumVarArray> &ys = my_y_vars[rId];
			vector<IloNumVarArray> &xs = my_x_vars[rId];
			for (int u=0; u<my_proc_num; u++) {
				IloNumVarArray y(env);
				IloNumVarArray x(env);
				for (int k=0; k<my_req_num; k++) {
					y.add(IloNumVar(env, 0, 1, ILOINT));
					x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
				}
				ys.push_back(y);
				xs.push_back(x);
			}
		}
		//		cout << "Stop populating my variables" << endl;

		/* Impose constraints */
		/* Constraint 3: for each request to a resource, sum 
		 * of corresponding y variables is at most 1
		 * Constraint 2 then follows
		 */
		//		cout << "Start imposing constraint 3" << endl;
		map<ResourceID, vector<IloNumVarArray> >::iterator y_it = my_y_vars.begin();
		for (; y_it != my_y_vars.end(); y_it++) {
			ResourceID rId = y_it->first;			
			vector<IloNumVarArray> &var_arrays = y_it->second;
			unsigned int req_num = myRes[rId]->requestNum;
			for (int k=0; k<req_num; k++) {
				IloExpr expr(env);
				for (int u=0; u<my_proc_num; u++) {
					expr += var_arrays[u][k];
				}
				cons.add(expr <= 1);
				expr.end(); /* Important: must free it */
			}
		} /* End constraint 3 */
		//		cout << "End imposing constraint 3" << endl;
		
		/* Constraint 5: products of mismatched y and x variables 
		 * are all zero
		 */
		//		cout << "Start imposing constraints 5 and 6" << endl;
		y_it = my_y_vars.begin();
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
		} /* End constraints 5 & 6 */
		//		cout << "Stop imposing constraints 5 and 6" << endl;

		/* Constraint 8 (FIFO): at most 1 request from each other
		 * processor can block me
		 */
		//		cout << "Start imposing constraint 8" << endl;
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

		y_it = my_y_vars.begin();
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
		//		cout << "Stop imposing constraint 8" << endl;
		
		/* Add constraints to Cplex model */
		model.add(cons);

		/* Objective function */
		//		cout << "Building objective function" << endl;
		IloExpr obj(env);
		inter_it = x_vars.begin();
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

		y_it = my_y_vars.begin();
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
		//		cout << "Finish building objective function" << endl;

		model.add(IloMaximize(env, obj));
		obj.end();
		IloCplex cplex(model);

		if (!cplex.solve()) {
			env.error() << "Failed to optimize" << endl;
			throw(-1);
		} else {
			cout << "SOLVE SUCCESSFULLY" << endl;
		}

		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value  = " << cplex.getObjValue() << endl;

		cplex.exportModel("fifo.lp");
	} catch (IloException &e) {
		cerr << "Ilog concert exception: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception: " << endl;
	}

	cout << endl;
	env.end();
}

void blocking_analysis(TaskSet* tset) {
	
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
