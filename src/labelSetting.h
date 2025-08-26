#ifndef LABELSETTING_H
#define LABELSETTING_H

#ifdef __cplusplus
extern "C" {
#endif
#include "cli_parse.h"
#include "instance.h"
#include "data_structs.h"
#include "item.h"
#include "ts_heur.h"
#include <string.h>
#ifdef __cplusplus
}
#endif

#define EX 1e-8
#define MAXMUMSOLTIME 3600

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <stack>
#include <bitset>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include <thread>
#include <functional>
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>
#include "gurobi_c++.h"

/*parameters*/
extern int g_tolItemNum;							//the total number of items
extern int g_firstDim;								//the first dimension of the bucket
extern int secondDim;								//the second dimension of the bucket
extern double g_bestLB;								//the best lower bound		
extern int g_secondDimDiv;							//the divisor for second dimension of the bucket
extern double g_smallest_bDiva;						//the smallest ratio of a/b for the present instance	
extern int g_dualBoundRoundNum;						//the roundNum when executing the dual bound calculation
extern int g_maximumTimeDualBound;					//the time limitation for executing the dual bound calculation
extern bool g_solveKnapsack;						//if solve the knapsack problem

extern Instance g_instance;							//g_instance data				

/*results*/
extern int g_tolGeneratedLabel;

/*bucket index*/
extern std::vector<int> g_item_to_firstIdx;			//map item to item bucket index
extern std::vector<int> g_weight_to_secondIdx;			//map weight to capacity bucket index
extern int g_domWidth;									//the width of buckets dominance check
extern bool g_pInteger;								//remark if parameter p is integer 
extern bool g_aInteger;								//remark if parameter a is integer 

/*error control*/
extern bool g_errorControl;							//if the error control is enabled


typedef struct {
	size_t index;		//index of items in the orginstances
	size_t pos;			//index of deleted items in the g_instances
	double sum_p_ab;
	double profit;
} ItemIndex;

//The label class
class MyLabel {
public:
	MyLabel() {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = 0;
		parentLab = nullptr;
		//bitSet.resize(g_tolWrd);
	}
	/*MyLabel(int sr3Size) {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = 0;
		parentLab = nullptr;
		sr3Vis.resize(sr3Size);
	}*/
	MyLabel(double& binNumDual) {
		lastItem = -1;
		tolWeight = 0;
		sum_a = 0;
		sum_b = 0;
		tolProfit = binNumDual;
		parentLab = nullptr;
		//bitSet.resize(g_tolWrd);
	}
	MyLabel(const MyLabel* preLab) {
		lastItem = preLab->lastItem;
		tolWeight = preLab->tolWeight;
		sum_a = preLab->sum_a;
		sum_b = preLab->sum_b;
		tolProfit = preLab->tolProfit;
		parentLab = preLab->parentLab;
		itemSet = preLab->itemSet;
		//sr3Vis = preLab->sr3Vis;
		//bitSet = preLab->bitSet;
	}

	//compare two vectors with elements only take 0 or 1
	bool areEqual(const std::vector<int>& a, const std::vector<int>& b) {
		std::bitset<1024> ba, bb;
		for (size_t i = 0; i < a.size(); ++i) {
			ba[i] = a[i];
			bb[i] = b[i];
		}
		return ba == bb;
	}

	bool operator==(const MyLabel& other) {
		return areEqual(this->itemSet, other.itemSet);
	}
public:
	int lastItem;           //the last inserted item in the bag
	double sum_a;			//the summation of a part
	double sum_b;			//the summation of b part
	double tolWeight;       //the total weight of the items in the bag
	double tolProfit;			//the total profit of the items in the bag
	MyLabel* parentLab;		//the parent label of the current label
	std::vector<int> itemSet;	//the inserted items in the bag, Size = g_tolItemNum, value = 1 if a item is inserted
	//std::vector<int> sr3Vis;  //the visiting time to each sr3
	//vector<int> bitSet;     //the bit projection of the inserted items in the bag
};
using Bin = MyLabel;		// rename MyLabel as Bin

//the best kanpsack solution
class KnapsackSol {
public:

public:
	MyLabel* bestLab;			//the best label
	std::string bestItemSet;	    //the best item set
};

/*data structure for the BP*/
//The stage information
class Stage {
public:

public:
	std::unordered_map<int, int> stageCandidates;	//the items included in branch set, the second point to the together set index
	std::vector<std::vector<int>> togetherSet;		//the corresponding toghther set with stageCandidates.
	std::unordered_map<int, std::vector<int>> separationPairs;//record the separation relationships
	std::vector<std::vector<int>> options;			//all the combinations of visit choice
	std::vector<std::vector<int>> SR3Dual;			//the SR3 duals corresponding to the options, <index>
};

//The branch pattern, the item should be packed in the bin with weight in range "branchRange"
class BranchPattern {
public:
	bool operator ==(const BranchPattern& b) {
		return branchItem == b.branchItem &&
			branchRange == b.branchRange;
	}
	BranchPattern() {
		branchRange = { 0, g_instance.capacity };
		existCols = 0;
	}
public:
	int branchItem;						//the item used to branch
	std::pair<double, double> branchRange;	//the capacity range of the bin the item is packed
	int existCols;						//the item exists in multiple bins (number)
};

//The branch arc
class BranchInfo {
public:
	bool operator ==(const BranchInfo& b) {
		return itemPair == b.itemPair;
	}
public:
	BranchPattern prePattern;			//the branch pattern
	std::pair<int, int> itemPair;		//the branch item pair, i < j
	bool same_diff;						//0 denotes the same branch, 1 is the different branch
};


class BranchNode {
public:
	BranchNode(
		BranchNode* parentNode
	) {
		lowerBound = parentNode->lowerBound;
		branchArcs = parentNode->branchArcs;
		branchPatterns = parentNode->branchPatterns;
		stages = parentNode->stages;
		allCols.resize(parentNode->allCols.size());
		for (int i = 0; i < parentNode->allCols.size(); ++i)
			allCols[i] = new Bin(parentNode->allCols[i]);
		//preSolution = pSol; //dont need to be inherit
		SR3s = parentNode->SR3s;
		preSolution.clear();
		nodeDepth = parentNode->nodeDepth + 1;
	}
	BranchNode() {
		lowerBound = 0;
		nodeDepth = 0;
	}
	//deconstruct prenode, release memory
	~BranchNode() {
		for (auto& e : allCols) {
			delete e;
			e = nullptr;
		}
		allCols.clear();
	}

public:
	/*basic information*/
	double lowerBound;								//the lower bound
	std::vector<Bin*> allCols;						//all the remaining columns

	/*branch information*/
	std::vector<BranchInfo> branchArcs;			//record of branch arcs
	std::unordered_map<int, BranchPattern> branchPatterns;//record of branch patterns
	std::vector<Stage> stages;					//stages derived from the branch arcs
	std::unordered_map<int, double> preSolution;//the solution obtained from last CG
	int nodeDepth;								//the depth of the node in the search tree
	/*cut information*/
	std::vector<std::unordered_set<int>> SR3s;			//the SR3 found for the present node
};

//solution
extern std::vector<Bin*> g_bestSol_binpacking;
extern double g_heuristicTime;					//the time spent on heuristic to solve submodular Bin packing
extern double g_BPTime;							//the total time for BP


//control the decimals to keep of the dual variables
extern double g_controlNum;
extern int g_maxItemNum;						//maximum number of items can be contained in one bin
extern int g_binpackingUB;						//the upper bound of the bin packing problem
extern int g_upper_lower_dual;					//0 indicates using upper bounds(don't truncating), 1 indicates using lower bounds(truncating dual variables), 2 indicates using original dual variables
extern int g_dual_bound_parameter;				//0 denote using fixed parameters and 1 denotes using variable parameters
extern int g_maxNumSR3Comb;						//the maximum number of SR3s to add originally 10000, set to 0 to mask it
extern int g_maxNumDIs;							//the maximum number of DIs to add, set to 0 means we don't use DIs


/*function decleration*/
//convert vector to string
std::string JoinVector(const std::vector<int>& vec);

// instance data������,�Ӵ�С����
void sort_instance_by_p_ab(Instance& preInstances, ItemIndex*& indices);

void sort_instance_by_weight(Instance& instance, ItemIndex*& indices);

//design the labelsetting algorithm to solve the submodular knapsack problem
int LabelSettingSolveKnapsack(Args& args);

//label extention
void LabelExtention(Instance& preInstances, MyLabel* parentLab, MyLabel* preLab, int item);

//completion bound
int CompletionBound(
	MyLabel* lab,
	int prefirstIdx,
	DblMatrix* ub_matr
);

//The whole dominance logic
bool DominanceLogic(
	int preSecondIdx,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>*& newExtended,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>*& oldExtended,
	MyLabel* preLab,
	bool newLabelFalg,
	bool heuDom);

void JgeDominance(
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>& newExtended,
	std::vector<std::multimap<double, MyLabel*, std::greater<double>>>& oldExtended,
	MyLabel* lab,
	std::vector<MyLabel*>& dominatedOldLabs,
	bool newLabelFlag,
	bool heuDom,
	bool dualCalFlag);

//truncate a and b to speed up the algorithm and get dual bound 
void CalculateDualbound(
	Instance& g_Instance,
	const std::vector<double>& weightRec,
	DblMatrix& ub_matr,
	std::vector<double> SR3Duals,
	std::vector<std::unordered_set<int>> SR3s,
	std::vector<double>& preciseDuals,
	const ItemIndex* indicesRec,
	double& heuLableTime
);

int my_lin_relax(const Instance* RESTRICT const inst,
	DblMatrix* RESTRICT const output,
	int binCapacity,
	const std::vector<double>& weightRec,
	bool dualFlag);

//use the label algorithm as the subproblem of the bin packing problem
bool LabelSettingSolveSubproblem(std::vector<Bin*>& newCols, std::ofstream& outPut, const std::vector<Bin*>& allCols, double& preTime, std::vector<double>& preciseDuals, double& rcLB, bool heuPricingFlag);

//use the label algorithm as the subproblem of the bin packing problem, used for child nodes
bool LabelSettingSolveSubproblem_Child(
	std::vector<Bin*>& newCols,
	std::ofstream& outPut,
	const std::vector<Bin*>& allCols,
	std::vector<Stage>& preStages,
	double& preTime,
	std::vector<double>& preciseDuals,
	std::vector<double>& SR3Duals,
	std::vector<std::unordered_set<int>>& SR3s,
	double& rcLB,
	std::unordered_map<int, BranchPattern>& branchPatterns,
	bool heuPricingFlag);

//call the branch and price to solve the submodular bin packing problem
void Branch_and_pricing(Args& args);

//copy instance
void CopyInstance(Instance& preInstance, Instance& orgInstance);

//call the heuristic to solve the submodular bin packing problem
int HeuristicForBinpacking(std::vector<Bin*>& bestSol);

void AdjustInstancesPerStages(
	std::vector<Stage>& preStages,
	ItemIndex*& indicesRec,
	std::unordered_map<int, ItemIndex>& removedIndicesRec);

int SolveCompactPricingModel(
	std::vector<Bin*>& newCols,
	std::ofstream& outPut,
	std::vector<Stage>& preStages,
	double& preTime,
	double& binNumDual,
	std::vector<double>& preciseDuals
);

// call the gurobi solver to solve the compact SMBP model
int SolveCompactSMBPModel(Args& args, double binLB, int binUB);

#endif
