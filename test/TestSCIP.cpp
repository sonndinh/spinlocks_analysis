#include <iostream>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

int main(int argc, char** argv) {
	char *file_name = argv[1];
	SCIP *scip = NULL;
	SCIP_CALL( SCIPcreate(&scip) );

	// load default plugins
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
	
	// disable output to stdout
	SCIP_CALL( SCIPsetMessagehdlr(scip, NULL) );

	SCIP_CALL( SCIPreadProb(scip, file_name, NULL) );
	
	std::cout << "Solving problem" << std::endl;
	SCIP_CALL( SCIPsolve(scip) );

	//	std::cout << "Print best solution" << std::endl;
	//	SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );

	// get best solution
	SCIP_Real obj_val = SCIPgetPrimalbound(scip);
	std::cout << "Best solution for the problem: " << std::endl;
	std::cout << "Objective value: " << obj_val << std::endl;
	
	SCIP_CALL( SCIPfree(&scip) );

	BMScheckEmptyMemory();
	
	return SCIP_OKAY;
}
