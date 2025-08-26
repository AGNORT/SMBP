#include<iostream>
#include "compactModel.h"
#include "labelSetting.h"

using namespace std;


int main(int argc, char** argv) {
	//read arguments
	Args args = { 0 };
	parse_args(&args, argc, argv);

	/*solve submodular knapsack problem*/
	//call labelsetting algorithm to solve the submodular knapsack problem
	/*cout << "*****breadth first search******" << endl;
	LabelSettingSolveKnapsack(args);*/
	////call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
	//SolveCompactKnapsackModel(args);

	///*Solve SMBP*/
	////Solve the SMBP problem with Gurobi solver
	//SolveCompactSMBPModel(args, 0, 100);

	//Solve the SMBP problem with BP
	Branch_and_pricing(args);


	return 0;
}
