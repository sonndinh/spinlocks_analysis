#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

/**
 * Test CPLEX with the following problem (which returns 3.3e-7
 * for objective value)
 *
Maximize
 obj: 0 x9 + 0 x10 + 0 x11 + 0 x12 + 0 x13 + 0 x14 + 0 x15 + 0 x16 + 11 x17
      + 11 x18 + 11 x19 + x20 + [ 22 x1 * x16 + 22 x2 * x12 + 22 x3 * x15
      + 22 x4 * x11 + 22 x5 * x14 + 22 x6 * x10 + 22 x7 * x13 + 22 x8 * x9
      ] / 2
Subject To
 c1: x1 + x2 <= 1
 c2: x3 + x4 <= 1
 c3: x5 + x6 <= 1
 c4: x7 + x8 <= 1
 c5: - 2 x1 - 2 x3 - 2 x5 - 2 x7 + x17 + x18 + x19 <= 0
 q1: [ x1 * x12 + x3 * x11 + x5 * x10 + x7 * x9 ] <= 0
 q2: [ x2 * x16 + x4 * x15 + x6 * x14 + x8 * x13 ] <= 0
 q3: [ x1 * x16 + x3 * x15 + x5 * x14 + x7 * x13 ] <= 0
 q4: - x1 - x3 - x5 - x7 + [ x2 * x12 + x4 * x11 + x6 * x10 + x8 * x9 ] <= 0
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
 0 <= x8 <= 1
 1e-8 <= x9 <= 1
 1e-8 <= x10 <= 1
 1e-8 <= x11 <= 1
 1e-8 <= x12 <= 1
 1e-8 <= x13 <= 1
 1e-8 <= x14 <= 1
 1e-8 <= x15 <= 1
 1e-8 <= x16 <= 1
 1e-8 <= x17 <= 1
 1e-8 <= x18 <= 1
 1e-8 <= x19 <= 1
      x20 = 0
Binaries
 x1  x2  x3  x4  x5  x6  x7  x8 
*/

static void populate(IloModel model, IloNumVarArray var, IloRangeArray con);

int main(int argc, char** argv) {
	IloEnv env;
	try {
		IloModel model(env);
		IloNumVarArray var(env);
		IloRangeArray con(env);

		populate(model, var, con);

		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize" << endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, var);
		env.out() << "Vals: " << vals << endl;
		
		cplex.exportModel("testqcp.lp");

	} catch (IloException &e) {
		cerr << "Concert exception caught: " << e << endl;
	} catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();
	return 0;
}

static void populate(IloModel model, IloNumVarArray x, IloRangeArray c) {
	IloEnv env = model.getEnv();

	/* Add 8 binary variables: x1 to x8 */
	for (int i=0; i<8; i++) {
		x.add(IloNumVar(env, 0, 1, ILOINT));
	}

	/* Add 11 real variables: x9 to x19 */
	for (int i=0; i<11; i++) {
		x.add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
		//x.add(IloNumVar(env, 1e-8, 1.0, ILOFLOAT));
	}

	/* Add x20 (not sure what it is used for) */
	//	x.add(IloNumVar(env, 0.0, 0.0));

	c.add(1 <= x[0] + x[1] <= 1);
	c.add(1 <= x[2] + x[3] <= 1);
	c.add(1 <= x[4] + x[5] <= 1);
	c.add(1 <= x[6] + x[7] <= 1);

	c.add( -2*x[0] - 2*x[2] - 2*x[4] - 2*x[6] + 
		   x[16] + x[17] + x[18] <= 0);

	/* q1 constraint: add (x*x - x) for each binary variable in an attempt to make 
	 * it positive semi-definite 
	 */
	//	c.add(x[0]*x[11] + x[2]*x[10] + x[4]*x[9] + x[6]*x[8] <= 0);
	//	x[0]*x[0] - x[0] + x[2]*x[2] - x[2] + x[4]*x[4] - x[4] + x[6]*x[6] - x[6] +
	//		x[1]*x[1] - x[1] + x[3]*x[3] - x[3] + x[5]*x[5] - x[5] + x[7]*x[7] - x[7] <= 0);

	/* q2 constraint : get rid of mismatched combination
	 */
	//	c.add(x[1]*x[15] + x[3]*x[14] + x[5]*x[13] + x[7]*x[12] <= 0);
		  //		  x[0]*x[0] - x[0] + x[2]*x[2] - x[2] + x[4]*x[4] - x[4] + x[6]*x[6] - x[6] +
		  //		  x[1]*x[1] - x[1] + x[3]*x[3] - x[3] + x[5]*x[5] - x[5] + x[7]*x[7] - x[7] <= 0);

	/* q3 constraint: get rid of mismatched combination
	 */
	//	c.add(x[0]*x[15] + x[2]*x[14] + x[4]*x[13] + x[6]*x[12] <= 0);
		  //		  x[0]*x[0] - x[0] + x[2]*x[2] - x[2] + x[4]*x[4] - x[4] + x[6]*x[6] - x[6] +
		  //		  x[1]*x[1] - x[1] + x[3]*x[3] - x[3] + x[5]*x[5] - x[5] + x[7]*x[7] - x[7] <= 0);

	/* q4 constraint: fifo constraint for contention inside the task
	 */
	//	c.add( -x[0] - x[2] - x[4] - x[6] + x[1]*x[11] + x[3]*x[10] +
	//		   x[5]*x[9] + x[7]*x[8] <=0);
	
	model.add(c);

	model.add(IloMaximize(env, 11*(x[16] + x[17] + x[18] + 
								   x[0]*x[15] + x[1]*x[11] + 
								   x[2]*x[14] + x[3]*x[10] +
								   x[4]*x[13] + x[5]*x[9] + 
								   x[6]*x[12] + x[7]*x[8] )));
								   //								   x[0]*x[0] - x[0] + x[1]*x[1] - x[1] + 
								   //								   x[2]*x[2] - x[2] + x[3]*x[3] - x[3] + 
								   //								   x[4]*x[4] - x[4] + x[5]*x[5] - x[5] + 
								   //								   x[6]*x[6] - x[6] + x[7]*x[7] - x[7] )));
}
