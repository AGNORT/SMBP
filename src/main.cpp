#include<iostream>
#include "compactModel.h"
#include "labelSetting.h"
#include <string.h>

using namespace std;


int main(int argc, char** argv) {
	/**********Read arguments*******/
	Args args = { 0 };
	parse_args(&args, argc, argv);


	/**********Solve SMKP*******/
	if (strcmp(args.problem_type, "SMKP") == 0) {
		cout << "Solve SMKP use ";
		if (strcmp(args.method, "MyMethod") == 0) {
			cout << "DP!" << endl;
			//call labelsetting algorithm to solve the submodular knapsack problem
			LabelSettingSolveKnapsack(args);
		}
		else{
			cout << "Gurobi!" << endl;
			//call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
			SolveCompactKnapsackModel(args);
		}
	}
	else { //args.problem_type == "SMBP"
		cout << "Solve SMBP use ";
		/*******Solve SMBP*********/
		if (strcmp(args.method, "MyMethod") == 0){
			cout << "BPC!" << endl;
			
			//set parameters for SMBP
			if (strcmp(args.algoCombination, "ULA-VD") == 0) {
				g_upper_lower_dual = 0;
				g_dual_bound_parameter = 1;
			}
			else if (strcmp(args.algoCombination, "ULA-FD") == 0) {
				g_upper_lower_dual = 0;
				g_dual_bound_parameter = 0;
			}
			else if (strcmp(args.algoCombination, "LLA-VD") == 0) {
				g_upper_lower_dual = 1;
				g_dual_bound_parameter = 0;
			}
			else {//"LLA-FD"
				g_upper_lower_dual = 1;
				g_dual_bound_parameter = 1;
			}
			cout << "Use algorithm combination: " << args.algoCombination << endl;
			//SR3 and DIs
			if (!args.SR3) {
				g_maxNumSR3Comb = 0;
				cout << "Not use SR3 inequalities!" << endl;
			}
			if (args.DIs) {
				g_maxNumDIs = 10000; //maximum add 10000 DIs
				cout << "Use DIs inequalities!" << endl;
			}

			//Solve the SMBP problem with BPC
			Branch_and_pricing(args);
		}
		else{
			cout << "Gurobi!" << endl;
			//Solve the SMBP problem with Gurobi solver
			SolveCompactSMBPModel(args, 0, 100);
		}
	}

	return 0;
}
