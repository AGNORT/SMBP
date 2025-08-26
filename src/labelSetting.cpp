/****************The second version****************/
#include "labelSetting.h"
#include "compactModel.h"

using namespace std;

Instance g_instance = { 0 };					//g_instance data

/*parameters*/
int g_tolItemNum = 0;							//the total number of items
int g_firstDim = 0;								//the first dimension of the bucket
int secondDim = 0;								//the second dimension of the bucket
double g_bestLB = 0;								//the best lower bound		
int g_secondDimDiv = 1;							//the divisor for second dimension of the bucket
double g_smallest_bDiva = 1e10;					//the smallest ratio of a/b for the present instance	

vector<double> g_profitRec;						//the profit record from each item

/*bucket index*/
std::vector<int> g_item_to_firstIdx;			//map item to item bucket index
std::vector<int> g_weight_to_secondIdx;			//map weight to capacity bucket index
int g_domWidth = 6;								//the width of buckets dominance check
bool g_pInteger = true;						//remark if parameter p is integer 
bool g_aInteger = true;						//remark if parameter a is integer 

/*results*/
int g_tolGeneratedLabel = 0;					//the total number of generated labels
int g_nonDominatedLabel = 0;					//the total number of non-dominated labels
int g_CBFathomLabel = 0;						//the number of labels fathomed by CB
int g_dominatedLabel = 0;						//the number of labels fathomed by CB

/*BP parameters*/
int g_maxAddColNumPerCG = 1000;					//the maximum added number of columns per CG iteration

/*error control*/
bool g_errorControl = false;					//if the error control is enabled

/*sort the instances according to the sum of a and b*/
// ascending sort
int compare_item_index_ab(const void* a, const void* b) {
	double diff = ((ItemIndex*)b)->sum_p_ab - ((ItemIndex*)a)->sum_p_ab;
	return (diff > 0) - (diff < 0);
}
int compare_item_index_prt(const void* a, const void* b) {
	double diff = ((ItemIndex*)b)->profit - ((ItemIndex*)a)->profit;
	return (diff > 0) - (diff < 0);
}
// descending sort
int compare_item_index_descending(const void* a, const void* b) {
	double diff = ((ItemIndex*)a)->sum_p_ab - ((ItemIndex*)b)->sum_p_ab;
	return (diff > 0) - (diff < 0);
}
// instance data������,�Ӵ�С����
void sort_instance_by_p_ab(Instance& preInstances, ItemIndex*& indices) {
	size_t n = preInstances.n_items;
	indices = (ItemIndex*)malloc(n * sizeof(ItemIndex));

	// ��ʼ����������
	for (size_t i = 0; i < n; ++i) {
		indices[i].index = i;
		indices[i].sum_p_ab = preInstances.p_ptr[i] / (preInstances.a_ptr[i] + sqrt(preInstances.b_ptr[i]));
		//indicesRec[i].sum_p_ab = 1 / (preInstances.a_ptr[i] + sqrt(preInstances.b_ptr[i]));
	}

	// ����
	qsort(indices, n, sizeof(ItemIndex), compare_item_index_ab);

	// �����µ����鲢�������������
	double* new_a = (double*)malloc(n * sizeof(double));
	double* new_b = (double*)malloc(n * sizeof(double));
	double* new_p = (double*)malloc(n * sizeof(double));
	double* new_pw = (double*)malloc(n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		size_t idx = indices[i].index;
		new_a[i] = preInstances.a_ptr[idx];
		new_b[i] = preInstances.b_ptr[idx];
		new_p[i] = preInstances.p_ptr[idx];
		new_pw[i] = preInstances.p_weight[idx];
	}

	// �滻ԭʼָ��
	free(preInstances.a_ptr);
	free(preInstances.b_ptr);
	free(preInstances.p_ptr);
	free(preInstances.p_weight);
	preInstances.a_ptr = new_a;
	preInstances.b_ptr = new_b;
	preInstances.p_ptr = new_p;
	preInstances.p_weight = new_pw;
}
// instance data������,����profit�Ӵ�С����
void sort_instance_by_profit(Instance& instance) {
	size_t n = instance.n_items;
	ItemIndex* indices = (ItemIndex*)malloc(n * sizeof(ItemIndex));

	// ��ʼ����������
	for (size_t i = 0; i < n; i++) {
		indices[i].index = i;
		indices[i].sum_p_ab = instance.p_ptr[i];
	}

	// ����
	qsort(indices, n, sizeof(ItemIndex), compare_item_index_prt);

	// �����µ����鲢�������������
	double* new_a = (double*)malloc(n * sizeof(double));
	double* new_b = (double*)malloc(n * sizeof(double));
	double* new_p = (double*)malloc(n * sizeof(double));
	double* new_pw = (double*)malloc(n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		size_t idx = indices[i].index;
		new_a[i] = instance.a_ptr[idx];
		new_b[i] = instance.b_ptr[idx];
		new_p[i] = instance.p_ptr[idx];
		new_pw[i] = instance.p_weight[idx];
	}

	// �滻ԭʼָ��
	instance.a_ptr = new_a;
	instance.b_ptr = new_b;
	instance.p_ptr = new_p;
	instance.p_weight = new_pw;

	free(indices);
}
// instance data sort according to a + rho*sqrt(b), ascending
void sort_instance_by_weight(Instance& instance, ItemIndex*& indices) {
	size_t n = instance.n_items;
	indices = (ItemIndex*)malloc(n * sizeof(ItemIndex));

	for (size_t i = 0; i < n; i++) {
		indices[i].index = i;
		indices[i].sum_p_ab = instance.a_ptr[i] + instance.rho * sqrt(instance.b_ptr[i]);
	}

	qsort(indices, n, sizeof(ItemIndex), compare_item_index_descending);

	double* new_a = (double*)malloc(n * sizeof(double));
	double* new_b = (double*)malloc(n * sizeof(double));
	double* new_p = (double*)malloc(n * sizeof(double));
	double* new_pw = (double*)malloc(n * sizeof(double));

	for (size_t i = 0; i < n; i++) {
		size_t idx = indices[i].index;
		new_a[i] = instance.a_ptr[idx];
		new_b[i] = instance.b_ptr[idx];
		new_p[i] = instance.p_ptr[idx];
		new_pw[i] = instance.p_weight[idx];
	}

	instance.a_ptr = new_a;
	instance.b_ptr = new_b;
	instance.p_ptr = new_p;
	instance.p_weight = new_pw;
}


//convert vector to string
std::string JoinVector(const std::vector<int>& vec) {
	std::ostringstream oss;
	for (size_t i = 0; i < vec.size(); ++i) {
		if (i != 0) oss << ",";
		oss << vec[i];
	}
	return oss.str();
}

////find suitable items to be inserted into the bag
//void FindSuitItems(MyLabel* preLab, vector<int>& itemToEx) {
//	double lastWeight = -1;
//	for (int j = preLab->lastItem + 1; j < g_tolItemNum; ++j) {
//		////Judge duplicated items
//		//if ((preLab->bitSet[g_wrd[j]] & g_bit[j]) == g_bit[j])
//		//	continue;
//
//		double preWeight = preLab->sum_a + g_instance.a_ptr[j] +
//			g_instance.rho * sqrt(preLab->sum_b + g_instance.b_ptr[j]);
//		if(
//			!itemToEx.empty() &&
//			j > preLab->lastItem + 1 &&
//			preWeight >= lastWeight &&
//			g_instance.p_ptr[j]  <= g_instance.p_ptr[j-1]
//			//g_instance.p_ptr[j] / preWeight <= (g_instance.p_ptr[j - 1] / lastWeight)
//		  )
//			break;
//
//		//judge capacity
//		if (preWeight <= g_instance.capacity) {
//			itemToEx.push_back(j);
//			lastWeight = preWeight;
//		}
//		else {//update itemSet and bitSet
//			if (j == preLab->lastItem + 1) {
//				preLab->lastItem = g_tolItemNum - 1;			//lastItem extension
//				//cout << "Hit lastItem extension" << endl;
//			}
//			//break;
//			//preLab->itemSet.push_back(j);
//			//preLab->bitSet[g_wrd[j]] |= g_bit[j];
//		}
//	}
//}

//call the gurobi solver to solve the linear relaxation of the knapsack problem
void GetExactLinearSol(Instance& instance, DblMatrix& ub_matr, vector<double>& weightRec) {
	//compute linear upper bound (dual bound)
	for (size_t i = 1; i <= g_instance.n_items; ++i) {
		for (size_t j = 1; j <= secondDim; ++j) {
			double dualBound = INFINITY;
			SolveLinearKnapsack(instance, i, weightRec[j], dualBound);
			dm_set(&ub_matr, i, j, INFINITY);
		}
	}
}

//label extention
void LabelExtention(
	Instance& preInstance,
	MyLabel* parentLab,
	MyLabel* preLab,
	int item
) {
	preLab->lastItem = item;
	preLab->sum_a += preInstance.a_ptr[item];
	preLab->sum_b += preInstance.b_ptr[item];
	preLab->tolWeight = preLab->sum_a + preInstance.rho * sqrt(preLab->sum_b);
	preLab->tolProfit += preInstance.p_ptr[item];
	preLab->parentLab = parentLab;
	////update SR3s, consider the dual variables
	//for (int i = 0; i < SR3s.size(); ++i) {
	//	if (SR3s[i].find(item) != SR3s[i].end()) {
	//		if (preLab->sr3Vis[i] == 1)
	//			preLab->tolProfit += SR3Duals[i];
	//		++preLab->sr3Vis[i];
	//		preLab->sr3Vis[i] %= 2;
	//	}
	//}
}

//Judge same label
bool JudgeEquality(MyLabel* lab1, MyLabel* lab2) {
	return abs(lab1->tolProfit - lab2->tolProfit) <= EX &&
		abs(lab1->tolWeight - lab2->tolWeight) <= EX &&
		abs(lab1->sum_a - lab2->sum_a) <= EX;
}

//dominance rules
bool JgeLabDominance(MyLabel* lab1, MyLabel* lab2, bool heuDom) {
	if (!heuDom) {//exact dominance rules
		if (!g_pInteger) {
			if (lab1->tolProfit < lab2->tolProfit - EX) return false;
		}
		else {
			if ((int)lab1->tolProfit < (int)lab2->tolProfit) return false;
		}
		//if (lab1->sum_b > lab2->sum_b + EX) return false; //old weaker version
		if (lab1->tolWeight > lab2->tolWeight + EX) return false;
		if (!g_aInteger) {
			if (lab1->sum_a > lab2->sum_a + EX) return false;
		}
		else {
			if ((int)lab1->sum_a > (int)lab2->sum_a) return false;
		}


		if (g_pInteger && !g_aInteger) {
			if (
				(int)lab1->tolProfit > (int)lab2->tolProfit ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				lab1->sum_a < lab2->sum_a - EX) return true;
		}
		else if (!g_pInteger && g_aInteger) {
			if (
				lab1->tolProfit > lab2->tolProfit + EX ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				(int)lab1->sum_a < (int)lab2->sum_a) return true;
		}
		else if (!g_pInteger && !g_aInteger) {
			if (
				lab1->tolProfit > lab2->tolProfit + EX ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				lab1->sum_a < lab2->sum_a - EX) return true;
		}
		else {//(g_pInteger && g_aInteger)
			if (
				(int)lab1->tolProfit > (int)lab2->tolProfit ||
				//lab1->sum_b < lab2->sum_b - EX || //old weaker version
				lab1->tolWeight < lab2->tolWeight - EX ||
				(int)lab1->sum_a < (int)lab2->sum_a) return true;
		}

		return false;
	}
	else {//heurisitc dominance rules
		if (lab1->tolProfit < lab2->tolProfit - EX) return false;
		if (lab1->tolWeight > lab2->tolWeight + EX) return false;

		if (
			lab1->tolProfit > lab2->tolProfit + EX ||
			lab1->tolWeight < lab2->tolWeight - EX
			) return true;
		return false;
	}
}

//completion bound
int CompletionBound(
	MyLabel* lab,
	int prefirstIdx,
	DblMatrix* ub_matr
) {
	////original completion bound
	//if (lab->tolProfit + g_profitRec[lab->lastItem] <= g_bestLB - EX)
	//	return 2;
	//return false;

	//strengthened completion bound
	//size_t i = ub_matr->cols - (lab->lastItem + 1);
	size_t i = ub_matr->cols - prefirstIdx - 1;
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];
	size_t j = ub_matr->rows - secondIdx - 1;

	if (i < g_firstDim) {
		double leftVal = dm_get(ub_matr, i, j) + lab->tolProfit;
		if (leftVal <= g_bestLB + EX) {//
			return true;
		}
	}
	return false;
	//size_t lo_j = ub_matr->rows - secondIdx - 1;
	//size_t up_j = ub_matr->rows - secondIdx - 1;
	//double best_compl = -1;
	//for (size_t j = lo_j; j <= up_j; ++j) {
	//	double fill = dm_get(ub_matr, i, j);
	//	double leftVal = fill + lab->tolProfit;
	//	if (best_compl < leftVal) {
	//		best_compl = leftVal;
	//	}
	//	if (leftVal > g_bestLB + EX) {
	//		return false;
	//	}
	//}
	////original completion bound
	//if (best_compl <= g_bestLB - EX)
	//	return true;
	//return false;
}

//label dominance
void JgeDominance(
	vector<vector<multimap<double, MyLabel*, greater<double>>>>& nonExtended,
	//vector<vector<multimap<double, MyLabel*, greater<double>>>>& extended,
	set<int>& firstNonBucketIdx,
	vector<set<int>>& secondNonBucketIdx,
	MyLabel* lab,
	int& labelCnt) {
	//find the dimension of the bucket
	int firstIdx = g_item_to_firstIdx[lab->lastItem];
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];

	////judge dominance
	bool dominanceFlag = false;
	auto idxIte1 = firstNonBucketIdx.begin();
	////for (int j = 0; j <= secondIdx; ++j) {
	//while (idxIte1 != firstNonBucketIdx.end()) {
	//	int i = *idxIte1;
	//	if (i > firstIdx) break;
	//	auto idxIte2 = secondNonBucketIdx[i].begin();
	//	while (idxIte2 != secondNonBucketIdx[i].end()) {
	//		int j = *idxIte2;
	//		if (j > secondIdx) break;
	//		/*if (!nonExtended[i][j].empty() &&
	//			JgeLabDominance(nonExtended[i][j].begin()->second, lab)) {
	//			dominanceFlag = true;
	//			delete lab;
	//			break;
	//		}
	//		if (!extended[i][j].empty() &&
	//			JgeLabDominance(extended[i][j].begin()->second, lab)) {
	//			dominanceFlag = true;
	//			delete lab;
	//			break;
	//		}*/
	//		auto domIite = nonExtended[i][j].begin();
	//		while (domIite != nonExtended[i][j].end()) {
	//			if (domIite->first < lab->tolProfit) break;
	//			if (JgeLabDominance(domIite->second, lab)) {
	//				dominanceFlag = true;
	//				delete lab;
	//				break;
	//			}
	//			++domIite;
	//		}
	//		if (dominanceFlag)break;
	//		//auto ite1 = extended[i][j].begin();
	//		//while (ite1 != extended[i][j].end()) {
	//		//	if (ite1->first < lab->tolProfit) break;
	//		//	if (JgeLabDominance(ite1->second, lab)) {
	//		//		dominanceFlag = true;
	//		//		delete lab;
	//		//		break;
	//		//	}
	//		//	++ite1;
	//		//}
	//		//if (dominanceFlag)break;
	//		++idxIte2;
	//	}
	//	if (dominanceFlag)break;
	//	++idxIte1;
	//}
	//if (!dominanceFlag) {
	//	for (int j = 0; j <= secondIdx; ++j) {
	//		auto ite1 = extended[j].begin();
	//		while (ite1 != extended[j].end()) {
	//			if (ite1->first < lab->tolProfit) break;
	//			if (JgeLabDominance(ite1->second, lab)) {
	//				dominanceFlag = true;
	//				delete lab;
	//				break;
	//			}
	//			++ite1;
	//		}
	//		if (dominanceFlag)break;
	//	}
	//}

	if (!dominanceFlag) {
		//use the present label to dominate other labels in the bucket
		idxIte1 = firstNonBucketIdx.find(firstIdx);
		while (idxIte1 != firstNonBucketIdx.end()) {
			int i = *idxIte1;
			auto idxIte2 = secondNonBucketIdx[i].find(secondIdx);
			auto tIte = idxIte2;
			while (idxIte2 != secondNonBucketIdx[i].end()) {
				int j = *idxIte2;
				if (!nonExtended[i][j].empty() &&
					lab->tolProfit < (--nonExtended[i][j].end())->first) {
					++idxIte2;
					continue;
				}
				/*if (idxIte2 != tIte) {
					if (!nonExtended[i][j].empty() &&
						nonExtended[i][j].begin()->first > 0 &&
						lab->tolProfit >= nonExtended[i][j].begin()->first) {
						auto domIite = nonExtended[i][j].begin();
						while (domIite != nonExtended[i][j].end()) {
							delete domIite->second;
							domIite = nonExtended[i][j].erase(domIite);
							--labelCnt;
						}
					}
					if (!extended[i][j].empty() &&
						lab->tolProfit >= extended[i][j].begin()->first) {
						auto domIite = extended[i][j].begin();
						while (domIite != extended[i][j].end()) {
							delete domIite->second;
							domIite = extended[i][j].erase(domIite);
							continue;
						}
					}
				}
				else {
					auto domIite = nonExtended[i][j].begin();
					while (domIite != nonExtended[i][j].end()) {
						if (JgeLabDominance(lab, domIite->second)) {
							delete domIite->second;
							domIite = nonExtended[i][j].erase(domIite);
							--labelCnt;
							continue;
						}
						++domIite;
					}
				}*/

				auto domIite = nonExtended[i][j].begin();
				while (domIite != nonExtended[i][j].end()) {
					if (JgeLabDominance(lab, domIite->second, true)) {
						delete domIite->second;
						domIite = nonExtended[i][j].erase(domIite);
						--labelCnt;
						//cout << "nonExtended dominanted" << endl;
						continue;
					}
					++domIite;
				}
				/*domIite = extended[j].begin();
				if (!extended[j].empty() &&
					lab->tolProfit < (--extended[j].end())->first) {
					++idxIte2;
					continue;
				}
				while (domIite != extended[j].end()) {
					if (JgeLabDominance(lab, domIite->second)) {
						delete domIite->second;
						domIite = extended[j].erase(domIite);
						continue;
					}
					++domIite;
				}*/
				++idxIte2;
			}
			++idxIte1;
		}

		//nonExtended[firstIdx][secondIdx].insert(make_pair(lab->tolProfit, lab));
		////update the non-empty bucket index
		//firstNonBucketIdx.insert(firstIdx);
		//secondNonBucketIdx[firstIdx].insert(secondIdx);
		//++labelCnt;
	}
}

//label dominance
void JgeDominance(
	vector<multimap<double, MyLabel*, greater<double>>>& newExtended,
	vector<multimap<double, MyLabel*, greater<double>>>& oldExtended,
	MyLabel* lab,
	vector<MyLabel*>& dominatedOldLabs,
	bool newLabelFlag = true,
	bool heuDom = false,
	bool dualCalFlag = false) {
	//find the dimension of the bucket
	int secondIdx = g_weight_to_secondIdx[(size_t)ceil(lab->tolWeight)];

	//use the present label to dominate other labels in the bucket
	int j = secondIdx;
	while (j < secondDim) {
		//dominate the labels in the new bucket
		if (!newExtended[j].empty() &&
			lab->tolProfit >= (--newExtended[j].end())->first) {
			auto ite = newExtended[j].begin();
			while (ite != newExtended[j].end()) {
				if (JgeLabDominance(lab, ite->second, heuDom)) {
					if (dualCalFlag) {
						delete ite->second;
						ite->second = nullptr;
					}
					else {
						if (j == secondIdx)
							dominatedOldLabs.push_back(ite->second);
						else {
							delete ite->second;
							ite->second = nullptr;
							++g_dominatedLabel;
						}
					}
					ite = newExtended[j].erase(ite);
					continue;
				}
				else {//Judge label equality
					if (JudgeEquality(lab, ite->second)) {
						if (dualCalFlag) {
							delete ite->second;
							ite->second = nullptr;
						}
						else {
							if (j == secondIdx)
								dominatedOldLabs.push_back(ite->second);
							else {
								delete ite->second;
								ite->second = nullptr;
								++g_dominatedLabel;
							}
						}
						ite = newExtended[j].erase(ite);
						continue;
					}
				}
				++ite;
			}
		}
		////new generated label cannot dominate the labels in the old bucket
		//if (newLabelFlag) {
		//	++idxIte2;
		//	continue;
		//}
		////dominate the labels in the old bucket
		//if (!oldExtended[j].empty() &&
		//	lab->tolProfit >= (--oldExtended[j].end())->first) {
		//	auto ite = oldExtended[j].begin();
		//	while (ite != oldExtended[j].end()) {
		//		if (JgeLabDominance(lab, ite->second)) {
		//			delete ite->second;
		//			ite->second = nullptr;
		//			ite = oldExtended[j].erase(ite);
		//			continue;
		//		}
		//		else {//Judge label equality
		//			if (JudgeEquality(lab, ite->second)) {
		//				delete ite->second;
		//				ite->second = nullptr;
		//				ite = oldExtended[j].erase(ite);
		//				continue;
		//			}
		//		}
		//		++ite;
		//	}
		//}

		++j;
		if (j - secondIdx > g_domWidth)
			break;
	}
}


//get the initial lower bound
double GetInitialLowerBound(double& bestLB) {
	Instance instance = g_instance;
	//��instance����profit�Ӵ�С����
	sort_instance_by_profit(instance);

	double sum_a = 0;
	double sum_b = 0;
	for (int i = 0; i < instance.n_items; ++i) {
		sum_a += instance.a_ptr[i];
		sum_b += instance.b_ptr[i];
		if (sum_a + instance.rho * sqrt(sum_b) > instance.capacity)
			break;
		else
			bestLB += instance.p_ptr[i];
	}

	//test
	//bestLB = 20900;

	MyLabel* preLab = new MyLabel();
	preLab->tolProfit = bestLB;
	return bestLB;
}

//initialize g_profitRec
void InitProfitRec() {
	vector<int> tmp(g_instance.n_items, 0);
	g_profitRec.resize(g_instance.n_items, 0);
	double tolProfit = 0;
	tmp[0] = g_instance.p_ptr[0];
	for (int i = 0; i < g_instance.n_items; ++i) {
		tolProfit += g_instance.p_ptr[i];
		if (i >= 1)
			tmp[i] = tmp[i - 1] + g_instance.p_ptr[i];
	}
	for (int i = 0; i < g_instance.n_items; ++i)
		g_profitRec[i] = tolProfit - tmp[i];
}

//The whole dominance logic
bool DominanceLogic(
	int preSecondIdx,
	vector<multimap<double, MyLabel*, greater<double>>>*& newExtended,
	vector<multimap<double, MyLabel*, greater<double>>>*& oldExtended,
	MyLabel* preLab,
	bool newLabelFalg = true,
	bool heuDom = false) {
	bool dominanceFlag = false;
	int j = preSecondIdx;
	while (j >= 0) {
		auto& preBkt = (*newExtended)[j];
		//use newExtended labels to dominante
		auto domIite = preBkt.begin();
		int cnt = 0;
		int preSize = preBkt.size();
		while (cnt < preSize) {
			if (domIite->first < preLab->tolProfit) break;
			//judge dominance
			if (JgeLabDominance(domIite->second, preLab, heuDom)) {
				dominanceFlag = true;
				break;
			}
			else {
				//judge equality
				if (JudgeEquality(domIite->second, preLab)) {
					dominanceFlag = true;
					//cout << "Equal label 1" << endl;
					break;
				}
			}
			++domIite;
			++cnt;
		}
		if (dominanceFlag) break;

		//old labels don't need to be checked from old labels, since it's checked during the geration of the new label
		if (newLabelFalg) {
			auto& preBkt = (*oldExtended)[j];
			auto domIite = preBkt.begin();
			int cnt = 0;
			int preSize = preBkt.size();
			//while (domIite != preBkt.end()) {
			while (cnt < preSize) {
				if (domIite->first < preLab->tolProfit) break;
				//judge dominance
				if (JgeLabDominance(domIite->second, preLab, heuDom)) {
					dominanceFlag = true;
					break;
				}
				else {
					//judge equality
					if (JudgeEquality(domIite->second, preLab)) {
						dominanceFlag = true;
						//cout << "Equal label 2" << endl;
						break;
					}
				}
				++domIite;
				++cnt;
			}
			if (dominanceFlag)break;
		}

		--j;
		if (preSecondIdx - j > g_domWidth) break;
	}

	return dominanceFlag;
}

//extend the non-dominanted label for all the operations in all the stages
void ExtendNDlabsToStages(
	vector<multimap<double, MyLabel*, greater<double>>>*& NDLabs,
	vector<MyLabel*>& dominatedOldLabs,
	unordered_map<int, ItemIndex>& removedIndicesRec,
	vector<Stage>& preStages,
	vector<double>& SR3Duals
) {
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	for (int stage = 0; stage < preStages.size(); ++stage) {
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*NDLabs)[sndBktIdx].empty()) continue;
			auto ite = (*NDLabs)[sndBktIdx].begin();
			while (ite != (*NDLabs)[sndBktIdx].end()) {
				if (ite->second->itemSet.empty())
					ite->second->itemSet.resize(g_instance.n_items + removedIndicesRec.size(), 0);
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				++g_tolGeneratedLabel;

				for (int t = 0; t < preStages[stage].options.size(); ++t) {
					unordered_set<int> preOperetion(
						preStages[stage].options[t].begin(),
						preStages[stage].options[t].end());
					/*for (auto& e : preStages[stage].options[t]) {
						int togetherIdx = preStages[stage].stageCandidates[e];
						if (togetherIdx != -1 &&
							!preStages[stage].togetherSet[togetherIdx].empty()) {
							preOperetion.insert(preStages[stage].togetherSet[togetherIdx].begin(),
								preStages[stage].togetherSet[togetherIdx].end() );
						}
					}*/
					//judge weight constraints
					double sum_a = ite->second->sum_a;
					double sum_b = ite->second->sum_b;
					double sum_p = ite->second->tolProfit;
					for (auto& e : preOperetion) {
						int itemIdx = removedIndicesRec[e].pos;
						sum_a += g_instance.a_ptr[itemIdx];
						sum_b += g_instance.b_ptr[itemIdx];
						sum_p += g_instance.p_ptr[itemIdx];
					}
					//considered dual variables of the SR3s
					for (auto& sr3DualIdx : preStages[stage].SR3Dual[t])
						sum_p += SR3Duals[sr3DualIdx];

					double preWeight = sum_a + g_instance.rho * sqrt(sum_b);
					if (preWeight <= g_instance.capacity + EX) {//judge capacity
						//label extension
						MyLabel* tmpLab = new MyLabel(ite->second);
						tmpLab->sum_a = sum_a;
						tmpLab->sum_b = sum_b;
						tmpLab->tolWeight = preWeight;
						tmpLab->tolProfit = sum_p;
						tmpLab->parentLab = ite->second;
						for (auto& e : preOperetion) {
							tmpLab->itemSet[e] = 1;		//record visit item information here
							tmpLab->lastItem = removedIndicesRec[e].pos;		//only avoid overflow
						}
						++g_tolGeneratedLabel;

						////test
						//if (ite->first <=1+EX && tmpLab->tolProfit > 1 + EX) {
						//	int a = 0;
						//	cout << "Hit here" << endl;
						//}

						//dominance check
						int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
						bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, NDLabs, tmpLab, true, true);

						if (dominanceFlag) {
							delete tmpLab;
							tmpLab = nullptr;
							++g_dominatedLabel;
						}
						else {
							//dominate other labels
							JgeDominance((*newExtended), (*NDLabs), tmpLab, dominatedOldLabs, true, true);
							//add the non dominated label into the new extended bucket
							(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
							//g_bestLB = max(g_bestLB, tmpLab->tolProfit);
						}
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*NDLabs)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, NDLabs, oldLab, false, true);
					if (dominanceFlag) {
						/*delete oldLab;
						oldLab = nullptr;*/
						dominatedOldLabs.push_back(oldLab);
						continue;
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
			}
		}

		delete NDLabs;
		NDLabs = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	}

	delete newExtended;
}


//use the heuristic labeling algorithm with heuristic dominance rules to solve the submodular knapsack problem
void LabelSettingHeuristic(
	DblMatrix& ub_matr,
	double& heuLableTime,
	const ItemIndex* indicesRec,
	vector<Stage>& preStages,
	vector<double>& SR3Duals,
	vector<unordered_set<int>>& SR3s,
	unordered_map<int, ItemIndex>& removedIndicesRec,
	multimap<double, KnapsackSol, greater<double>>& finalSols
	//ofstream& outPut
) {
	auto startTime = chrono::high_resolution_clock::now();

	//initialize the bucket
	vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
	dominatedOldLabs.reserve(1e+6);
	vector<int> nonDominatedLabsRec(g_firstDim);

	MyLabel* initLab = new MyLabel();
	oldExtended->begin()->insert({ 0, initLab });
	++g_tolGeneratedLabel;

	//label extention and dominance
	int currItem = 0;				//record the current item to extend
	int thisNonDominatedLabel = 0;
	for (int stage = 0; stage < g_firstDim; ++stage) {
		thisNonDominatedLabel = 0;
		/*bool checkedFlag = false;
		double p_low_iw = 0;
		int p_low_sndIdx = -1;*/
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*oldExtended)[sndBktIdx].empty()) continue;
			auto ite = (*oldExtended)[sndBktIdx].begin();
			while (ite != (*oldExtended)[sndBktIdx].end()) {
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
					dominatedOldLabs.push_back(ite->second);
					keepOld = false;
				}
				else {
					++g_tolGeneratedLabel;
				}

				double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
					g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
				if (preWeight <= g_instance.capacity + EX) {//judge capacity
					//label extension
					MyLabel* tmpLab = new MyLabel(ite->second);
					LabelExtention(g_instance, ite->second, tmpLab, currItem);
					++g_tolGeneratedLabel;

					//completion bound to fathom label
					if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
						++g_CBFathomLabel;
						delete tmpLab; tmpLab = nullptr;
						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, true);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
								continue;
							}
							(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
						}
						continue;
					}

					//dominance check
					int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
					bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab, true, true);

					if (dominanceFlag) {
						delete tmpLab;
						tmpLab = nullptr;
						++g_dominatedLabel;
					}
					else {
						//dominate other labels
						JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs, true, true);
						//add the non dominated label into the new extended bucket
						(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*oldExtended)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, true);
					if (dominanceFlag) {
						/*delete oldLab;
						oldLab = nullptr;*/
						dominatedOldLabs.push_back(oldLab);
						continue;
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
			}
			thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
			//use the present best dual variable to update
			if (!(*newExtended)[sndBktIdx].empty()) {
				if ((*newExtended)[sndBktIdx].begin()->first > g_bestLB) {
					KnapsackSol preSol;
					preSol.bestLab = new MyLabel((*newExtended)[sndBktIdx].begin()->second);
					//get the item set
					if (preSol.bestLab->itemSet.empty())
						preSol.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
					auto tmpLab = preSol.bestLab;
					while (tmpLab->parentLab != nullptr) {
						preSol.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
						tmpLab = tmpLab->parentLab;
					}
					finalSols.insert({ (*newExtended)[sndBktIdx].begin()->first, preSol });
					g_bestLB = max(g_bestLB, (*newExtended)[sndBktIdx].begin()->first);
				}
			}

			////use the worst labels in the newExtended to do the variable fixing
			//if (!(*newExtended)[sndBktIdx].empty() &&
			//	sndBktIdx >= 1 &&
			//	!checkedFlag) {
			//	MyLabel* thisLab = (--(*newExtended)[sndBktIdx].end())->second;
			//	p_low_iw = thisLab->tolProfit;
			//	p_low_sndIdx = sndBktIdx;
			//	checkedFlag = true;
			//}
		}
		////check variable fixing
		//if (checkedFlag && stage + 1 < g_firstDim) {
		//	if (p_low_iw + g_instance.p_ptr[stage + 1] +
		//		dm_get(&ub_matr, stage + 1, p_low_sndIdx) <= g_bestLB) {
		//		cout << "Hit variable fixing " << stage + 1 << endl;
		//	}
		//}

		++currItem;
		delete oldExtended;
		oldExtended = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		nonDominatedLabsRec[stage] = thisNonDominatedLabel;
	}

	//if Stages is not empty, do the final extensions
	if (!preStages.empty()) {//add new valid labels into oldExtended
		ExtendNDlabsToStages(oldExtended, dominatedOldLabs, removedIndicesRec, preStages, SR3Duals);
	}

	//record the best solution
	auto ite = --(*oldExtended).end();
	MyLabel* bestLab = nullptr;
	while (true) {
		if (!ite->empty()) {
			if (bestLab == nullptr)
				bestLab = ite->begin()->second;

			auto subIte = ite->begin();
			while (subIte != ite->end()) {
				/*if (1 - subIte->first >= 0)
					break;*/
				KnapsackSol preSol;
				preSol.bestLab = subIte->second;
				finalSols.insert({ subIte->second->tolProfit, preSol });
				subIte = ite->erase(subIte);
			}
		}
		if (ite == (*oldExtended).begin())
			break;
		--ite;
	}
	if (bestLab != nullptr && bestLab->lastItem > 0) {
		//get the item set
		vector<int> bestIS;
		MyLabel* tmpLab = bestLab;
		while (tmpLab->parentLab != nullptr) {
			bestIS.push_back(indicesRec[tmpLab->lastItem].index);
			tmpLab = tmpLab->parentLab;
		}
		sort(bestIS.begin(), bestIS.end());
		string bestItemSet = JoinVector(bestIS);
		finalSols.begin()->second.bestItemSet = bestItemSet;
		//get the item set for all the solutions
		for (auto& e : finalSols) {
			if (e.second.bestLab->itemSet.empty())
				e.second.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
			tmpLab = e.second.bestLab;
			while (tmpLab->parentLab != nullptr) {
				e.second.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
				tmpLab = tmpLab->parentLab;
			}
		}

		//store the best solution
		//KnapsackSol bestSol;
		//bestSol.bestItemSet = bestItemSet;
		//bestSol.bestLab = bestLab;
		//finalSols.insert({ bestLab->tolProfit, bestSol });
		if (preStages.empty())
			g_bestLB = max(g_bestLB, finalSols.begin()->first);
	}
	//free space
	for (auto& t : (*oldExtended)) {
		for (auto& e : t)
			delete e.second;
	}
	delete oldExtended;
	delete newExtended;
	g_dominatedLabel += dominatedOldLabs.size();
	for (auto& t : dominatedOldLabs)
		delete t;
	g_nonDominatedLabel = thisNonDominatedLabel;
	//cout << "The trend of non dominated labels in heuristic label algorithm: " << endl;
	//outPut << "The trend of non dominated labels in heuristic label algorithm: " << endl;
	//for (auto& e : nonDominatedLabsRec) {
	//	cout << e << "\t";
	//	outPut << e << "\t";
	//}
	cout << endl;
	//outPut << endl;
	auto endTime = chrono::high_resolution_clock::now();
	heuLableTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
}

//test solution
bool TestSolution(Instance& preInstance, vector<int>& sol_items) {
	double sum_a = 0;
	double sum_b = 0;
	double profit = 0;
	for (auto& e : sol_items) {
		if (e > 0) {
			sum_a += preInstance.a_ptr[e];
			sum_b += preInstance.b_ptr[e];
			profit += preInstance.p_ptr[e];
		}
	}
	if (sum_a + preInstance.rho * sqrt(sum_b) <= preInstance.capacity) {
		cout << "Solution feasible, with profit: " << profit << endl;
		return true;
	}
	else {
		cout << "Solution infeasible! " << endl;
		return false;
	}
}

//design the labelsetting algorithm to solve the submodular knapsack problem
int LabelSettingSolveKnapsack(Args& args) {
	auto startTime = chrono::high_resolution_clock::now();
	/*read g_instance*/
	g_instance = { 0 };
	if (instance_parse(&g_instance, args.input_file)) {
		cerr << "Error reading g_instance file" << endl;
		return 0;
	}
	if (g_instance.rho < 1 && !args.direct_rho)
		g_instance.rho = sqrt(g_instance.rho) / sqrt(1 - g_instance.rho);
	ofstream outPut(args.output_file, ios::app);
	//judge if p and a is integer or not
	for (int i = 0; i < g_instance.n_items; ++i) {
		if (abs(floor(g_instance.a_ptr[i]) - g_instance.a_ptr[i]) > EX) {
			g_aInteger = false;
		}
		if (abs(floor(g_instance.p_ptr[i]) - g_instance.p_ptr[i]) > EX) {
			g_pInteger = false;
		}
		double preRatio = g_instance.b_ptr[i] * 1.0 / (g_instance.a_ptr[i] * 1.0);
		if (preRatio < g_smallest_bDiva)
			g_smallest_bDiva = preRatio;
	}

	//specify the dimension of the bucket
	g_tolItemNum = g_instance.n_items;			//the first dimension of the bucket
	g_firstDim = g_instance.n_items / 1;		//the first dimension of the bucket
	secondDim = g_instance.capacity / g_secondDimDiv;       //the second dimension of the bucket

	/*use heuristic to get primal bound*/
	auto heuPrimal_startTime = chrono::high_resolution_clock::now();
	//get the initial lower bound
	multimap<double, KnapsackSol, greater<double>> finalSols;
	g_bestLB = 0;
	GetInitialLowerBound(g_bestLB);
	//sort the instance according to p/(a+sqrt(b))
	ItemIndex* indicesRec = nullptr;
	sort_instance_by_p_ab(g_instance, indicesRec);
	//TS heuristic
	vector<int> ts_sol_items(g_tolItemNum);
	double ts_sol = heuristic_solution(&g_instance, &ts_sol_items[0]);
	g_bestLB = max(g_bestLB, ts_sol);
	//TestSolution(g_instance, ts_sol_items);//test TS solution
	////initialize g_profitRec
	//InitProfitRec();
	double tsSol_sum_a = 0;
	double tsSol_sum_b = 0;
	vector<int> map_ts_sol_items(g_tolItemNum);
	for (int i = 0; i < g_tolItemNum; ++i) {
		if (ts_sol_items[i]) {
			tsSol_sum_a += g_instance.a_ptr[i];
			tsSol_sum_b += g_instance.b_ptr[i];
			map_ts_sol_items[indicesRec[i].index] = 1;
		}
	}
	MyLabel* tsSolLab = new MyLabel();
	tsSolLab->tolProfit = ts_sol;
	tsSolLab->itemSet = map_ts_sol_items;
	tsSolLab->tolWeight = tsSol_sum_a + g_instance.rho * sqrt(tsSol_sum_b);
	KnapsackSol preSol;
	preSol.bestLab = tsSolLab;
	finalSols.insert({ ts_sol, preSol });
	auto heuPrimal_endTime = chrono::high_resolution_clock::now();

	/*construct bucket index*/
	int itemGap = g_instance.n_items / g_firstDim;	//the item gap for each item bucket
	g_item_to_firstIdx.resize(g_instance.n_items + 1, 0);
	for (int j = 0; j <= g_firstDim; ++j) {
		int startItem = j * itemGap;
		int endItem = (j == g_firstDim - 1) ? g_instance.n_items : (j + 1) * itemGap - 1;
		endItem = min((size_t)endItem, g_instance.n_items);
		for (int w = startItem; w <= endItem; ++w)
			g_item_to_firstIdx[w] = j;
	}
	vector<double> weightRec(secondDim + 1, 0);
	int weightGap = g_instance.capacity / secondDim;	//the capacity gap for each capacity bucket
	g_weight_to_secondIdx.resize(g_instance.capacity + 1, 0);
	for (int j = 0; j <= secondDim; ++j) {
		int startWeight = j * weightGap;
		int endWeight = (j == secondDim - 1) ? g_instance.capacity : (j + 1) * weightGap - 1;
		endWeight = min(endWeight * 1.0, g_instance.capacity);
		weightRec[j] = endWeight;
		for (int w = startWeight; w <= endWeight; ++w)
			g_weight_to_secondIdx[w] = j;
	}

	/*compute linear upper bound (dual bound)*/
	auto LR_startTime = chrono::high_resolution_clock::now();
	DblMatrix ub_matr = { 0 };
	//dm_new(&ub_matr, g_instance.n_items + 1, g_instance.capacity + 1);
	//for (size_t i = 0; i < g_instance.n_items; ++i) {
	//	for (size_t j = 0; j < g_instance.capacity; ++j) {
	//		dm_set(&ub_matr, i, j, INFINITY);
	//	}
	//}
	dm_new(&ub_matr, g_instance.n_items + 1, secondDim + 1);
	for (size_t i = 0; i <= g_instance.n_items; ++i) {
		for (size_t j = 0; j <= secondDim; ++j) {
			dm_set(&ub_matr, i, j, INFINITY);
		}
	}
	my_lin_relax(&g_instance, &ub_matr, g_instance.capacity + 1, weightRec, false);
	//lin_relax(&g_instance, &ub_matr, g_instance.capacity + 1);
	//GetExactLinearSol(g_instance, ub_matr, weightRec);
	auto LR_endTime = chrono::high_resolution_clock::now();

	/*use the heuristic labeling algorithm to get better lower bound*/
	double heuLableTime = 0;
	unordered_map<int, ItemIndex> removedIndicesRec;
	vector<Stage> preStages;
	vector<double> sr3Duals;
	vector<unordered_set<int>> sr3s;
	vector<double> preciseDuals(g_instance.n_items, 0);
	for (int i = 0; i < g_instance.n_items; ++i)
		preciseDuals[i] = g_instance.p_ptr[i];

	double dualBoundTime = 0;
	/*CalculateDualbound(g_instance, weightRec, ub_matr, sr3Duals, sr3s, preciseDuals, indicesRec, dualBoundTime);
	cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;*/

	LabelSettingHeuristic(ub_matr, heuLableTime, indicesRec, preStages, sr3Duals, sr3s, removedIndicesRec, finalSols);
	cout << "heuristic labeling get lower bound: " << g_bestLB << endl;
	cout << "heuristic labeling use time: " << heuLableTime << "s" << endl;
	int HLA_sol = g_bestLB;

	//// truncating a and b to get the dual bound
	//double dualBoundTime = 0;
	/*CalculateDualbound(g_instance, weightRec, ub_matr, sr3Duals, sr3s, preciseDuals, indicesRec, dualBoundTime);
	cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;*/

	auto exactDP_startTime = chrono::high_resolution_clock::now();
	//initialize the bucket
	vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
	dominatedOldLabs.reserve(1e+6);
	vector<int> nonDominatedLabsRec(g_firstDim);

	MyLabel* initLab = new MyLabel();
	oldExtended->begin()->insert({ 0, initLab });
	++g_tolGeneratedLabel;

	//label extention and dominance
	int currItem = 0;				//record the current item to extend
	bool timeLimitFlag = false;
	int thisNonDominatedLabel = 0;
	for (int stage = 0; stage < g_firstDim; ++stage) {
		thisNonDominatedLabel = 0;
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*oldExtended)[sndBktIdx].empty()) continue;
			auto ite = (*oldExtended)[sndBktIdx].begin();
			while (ite != (*oldExtended)[sndBktIdx].end()) {
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
					dominatedOldLabs.push_back(ite->second);
					keepOld = false;
				}
				else {
					++g_tolGeneratedLabel;
				}

				double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
					g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
				if (preWeight <= g_instance.capacity + EX) {//judge capacity
					//label extension
					MyLabel* tmpLab = new MyLabel(ite->second);
					LabelExtention(g_instance, ite->second, tmpLab, currItem);
					++g_tolGeneratedLabel;

					//completion bound to fathom label
					if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
						++g_CBFathomLabel;
						delete tmpLab; tmpLab = nullptr;
						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
								continue;
							}

							(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
						}
						continue;
					}

					//dominance check
					int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
					bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab);

					if (dominanceFlag) {
						delete tmpLab;
						++g_dominatedLabel;
					}
					else {
						//dominate other labels
						JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs);
						//add the non dominated label into the new extended bucket
						(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*oldExtended)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
					if (dominanceFlag) {
						/*delete oldLab;
						oldLab = nullptr;*/
						dominatedOldLabs.push_back(oldLab);
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
				////test maximum solution time
				//auto endTime = chrono::high_resolution_clock::now();
				//if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
				//	cout << "The maximum solution time of my DP is reached!" << endl;
				//	outPut << "The maximum solution time of my DP is reached!" << endl;
				//	timeLimitFlag = true;
				//	break;
				//}
			}
			//if (timeLimitFlag) break;
			thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
			////use the present best dual variable to update
			//if (!(*newExtended)[sndBktIdx].empty()) {
			//	KnapsackSol preSol;
			//	preSol.bestLab = new MyLabel((*newExtended)[sndBktIdx].begin()->second);
			//	//get the item set
			//	if (preSol.bestLab->itemSet.empty())
			//		preSol.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
			//	auto tmpLab = preSol.bestLab;
			//	while (tmpLab->parentLab != nullptr) {
			//		preSol.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
			//		tmpLab = tmpLab->parentLab;
			//	}
			//	finalSols.insert({ (*newExtended)[sndBktIdx].begin()->first, preSol });
			//	g_bestLB = max(g_bestLB, (*newExtended)[sndBktIdx].begin()->first);
			//}
		}
		auto endTime = chrono::high_resolution_clock::now();
		if (g_BPTime + chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
			cout << "The maximum solution time of my DP is reached!" << endl;
			outPut << "The maximum solution time of my DP is reached!" << endl;
			timeLimitFlag = true;
		}
		if (timeLimitFlag) {
			auto ite = (*oldExtended).begin();
			for (auto& e : (*newExtended)) {
				for (auto& t : e)
					ite->insert({ t.first, t.second });
				++ite;
			}
			break;
		}
		++currItem;
		delete oldExtended;
		oldExtended = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		nonDominatedLabsRec[stage] = thisNonDominatedLabel;
	}

	//record the best solution
	auto ite = --(*oldExtended).end();
	MyLabel* bestLab = nullptr;
	while (true) {
		if (!ite->empty()) {
			bestLab = ite->begin()->second;
			ite->erase(ite->begin());
			break;
		}
		if (ite == (*oldExtended).begin())
			break;
		--ite;
	}
	string bestItemSet = finalSols.begin()->second.bestItemSet;
	if (bestLab != nullptr) {
		//get the item set
		vector<int> bestIS;
		MyLabel* tmpLab = bestLab;
		while (tmpLab->parentLab != nullptr) {
			bestIS.push_back(indicesRec[tmpLab->lastItem].index);
			tmpLab = tmpLab->parentLab;
		}
		sort(bestIS.begin(), bestIS.end());
		bestItemSet = JoinVector(bestIS);

		//store the best solution
		KnapsackSol bestSol;
		bestSol.bestItemSet = bestItemSet;
		bestSol.bestLab = bestLab;
		finalSols.insert({ bestLab->tolProfit, bestSol });
		g_bestLB = max(g_bestLB, finalSols.begin()->first);
	}
	if (bestItemSet.empty()) {
		for (int i = 0; i < finalSols.begin()->second.bestLab->itemSet.size(); ++i)
			if (finalSols.begin()->second.bestLab->itemSet[i])
				bestItemSet += to_string(i) + ",";
	}

	//free space
	for (auto& e : finalSols)
		delete e.second.bestLab;
	for (auto& t : (*oldExtended)) {
		for (auto& e : t)
			delete e.second;
	}
	delete oldExtended;
	delete newExtended;
	g_dominatedLabel += dominatedOldLabs.size();
	for (auto& lab : dominatedOldLabs)
		delete lab;
	free(indicesRec);
	free(g_instance.a_ptr);
	free(g_instance.b_ptr);
	free(g_instance.p_ptr);
	free(g_instance.p_weight);
	g_instance.a_ptr = g_instance.b_ptr = g_instance.p_ptr = g_instance.p_weight = nullptr;
	dm_free(&ub_matr);

	auto exactDP_endTime = chrono::high_resolution_clock::now();
	auto endTime = chrono::high_resolution_clock::now();
	double myLableTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
	cout << "The best solution obtained by my labelsetting algorithm is: " << g_bestLB << endl;
	cout << "The solution obtained by TS heuristic is: " << ts_sol << endl;
	cout << "The solution obtained by heuristic labeling algorithm is: " << HLA_sol << endl;
	cout << "The time for my labelsetting algorithm is: " << myLableTime << "s" << endl;
	cout << "The time for heuristic primal solution is: " << chrono::duration_cast<chrono::milliseconds>(heuPrimal_endTime - heuPrimal_startTime).count() / 1000.0 << "s" << endl;
	cout << "The time for linear relaxation (dual bound) is: " << chrono::duration_cast<chrono::milliseconds>(LR_endTime - LR_startTime).count() / 1000.0 << "s" << endl;
	cout << "The time for heuristic DP is: " << heuLableTime << "s" << endl;
	cout << "The time for exact dual bound is: " << dualBoundTime << "s" << endl;
	cout << "The time for final exact DP is: " << chrono::duration_cast<chrono::milliseconds>(exactDP_endTime - exactDP_startTime).count() / 1000.0 << "s" << endl;
	cout << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
	cout << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
	cout << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
	cout << "The number of dominated labels is: " << g_dominatedLabel << endl;
	cout << "The best item set is: " << bestItemSet << endl;
	//cout << "The trend of non dominated labels: " << endl;

	outPut << "\n\n****The resulte of my labelsetting algorithm****" << g_bestLB << endl;
	outPut << "The best solution obtained by my labelsetting algorithm is: " << g_bestLB << endl;
	outPut << "The solution obtained by TS heuristic is: " << ts_sol << endl;
	outPut << "The solution obtained by heuristic labeling algorithm is: " << HLA_sol << endl;
	outPut << "The time for my labelsetting algorithm is: " << myLableTime << "s" << endl;
	outPut << "The time for heuristic primal solution is: " << chrono::duration_cast<chrono::milliseconds>(heuPrimal_endTime - heuPrimal_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The time for linear relaxation (dual bound) is: " << chrono::duration_cast<chrono::milliseconds>(LR_endTime - LR_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The time for heuristic DP is: " << heuLableTime << "s" << endl;
	outPut << "The time for exact dual bound is: " << dualBoundTime << "s" << endl;
	outPut << "The time for final exact DP is: " << chrono::duration_cast<chrono::milliseconds>(exactDP_endTime - exactDP_startTime).count() / 1000.0 << "s" << endl;
	outPut << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
	outPut << "The best item set is: " << bestItemSet << endl;
	outPut << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
	outPut << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
	outPut << "The number of dominated labels is: " << g_dominatedLabel << endl;
	//outPut << "The trend of non dominated labels: " << endl;
	//for (auto& e : nonDominatedLabsRec) {
	//	cout << e << "\t";
	//	outPut << e << "\t";
	//}
	cout << endl;
	outPut << endl;
	outPut.close();

	return g_bestLB;
}

//test precise reduced cost
bool TestPreciseReducedCost(
	MyLabel* preLab,
	vector<double>& preciseDual,
	vector<unordered_set<int>>& SR3s,
	vector<double>& sr3Duals
) {
	double preProfit = 0;
	for (int i = 0; i < preLab->itemSet.size(); ++i) {
		if (preLab->itemSet[i])
			preProfit += preciseDual[i];
	}
	//consider the SR3 duals
	int cnt = 0;
	for (auto& sr3 : SR3s) {
		int node1 = *sr3.begin();
		int node2 = *(++sr3.begin());
		int node3 = *(++(++sr3.begin()));
		if (
			(preLab->itemSet[node1] && preLab->itemSet[node2]) ||
			(preLab->itemSet[node1] && preLab->itemSet[node3]) ||
			(preLab->itemSet[node2] && preLab->itemSet[node3])
			)
			preProfit += sr3Duals[cnt];
		++cnt;
	}
	return 1 - preProfit < 0 - 1e-6;
}

//use the label algorithm as the subproblem of the bin packing problem, used for root node only
bool LabelSettingSolveSubproblem(
	vector<Bin*>& newCols,
	ofstream& outPut,
	const vector<Bin*>& allCols,
	double& preTime,
	vector<double>& preciseDuals,
	double& rcLB,
	bool heuPricingFlag) {
	try {
		auto startTime = chrono::high_resolution_clock::now();
		//specify the dimension of the bucket
		g_firstDim = g_instance.n_items / 1;		//the first dimension of the bucket
		secondDim = g_instance.capacity / g_secondDimDiv;       //the second dimension of the bucket

		/*use heuristic to get primal bound*/
		//get the initial lower bound
		multimap<double, KnapsackSol, greater<double>> finalSols;
		g_bestLB = 0;
		GetInitialLowerBound(g_bestLB);
		//sort the instance according to p/(a+sqrt(b))
		ItemIndex* indicesRec = nullptr;
		sort_instance_by_p_ab(g_instance, indicesRec);
		//TS heuristic
		vector<int> ts_sol_items(g_tolItemNum);
		double ts_sol = heuristic_solution(&g_instance, &ts_sol_items[0]);
		double tsSol_sum_a = 0;
		double tsSol_sum_b = 0;
		vector<int> map_ts_sol_items(g_tolItemNum);
		for (int i = 0; i < g_tolItemNum; ++i) {
			if (ts_sol_items[i]) {
				tsSol_sum_a += g_instance.a_ptr[i];
				tsSol_sum_b += g_instance.b_ptr[i];
				map_ts_sol_items[indicesRec[i].index] = 1;
			}
		}
		MyLabel* tsSolLab = new MyLabel();
		tsSolLab->tolProfit = ts_sol;
		tsSolLab->itemSet = map_ts_sol_items;
		tsSolLab->tolWeight = tsSol_sum_a + g_instance.rho * sqrt(tsSol_sum_b);
		KnapsackSol preSol;
		preSol.bestLab = tsSolLab;
		finalSols.insert({ ts_sol, preSol });
		//the best primal bound is at least equal to 1
		g_bestLB = max(max(g_bestLB, 1.0), ts_sol);

		/*construct bucket index*/
		int itemGap = g_instance.n_items / g_firstDim;	//the item gap for each item bucket
		g_item_to_firstIdx.resize(g_instance.n_items + 1, 0);
		for (int j = 0; j <= g_firstDim; ++j) {
			int startItem = j * itemGap;
			int endItem = (j == g_firstDim - 1) ? g_instance.n_items : (j + 1) * itemGap - 1;
			endItem = min((size_t)endItem, g_instance.n_items);
			for (int w = startItem; w <= endItem; ++w)
				g_item_to_firstIdx[w] = j;
		}
		vector<double> weightRec(secondDim + 1, 0);
		int weightGap = g_instance.capacity / secondDim;	//the capacity gap for each capacity bucket
		g_weight_to_secondIdx.resize(g_instance.capacity + 1, 0);
		for (int j = 0; j <= secondDim; ++j) {
			int startWeight = j * weightGap;
			int endWeight = (j == secondDim - 1) ? g_instance.capacity : (j + 1) * weightGap - 1;
			endWeight = min(endWeight * 1.0, g_instance.capacity);
			weightRec[j] = endWeight;
			for (int w = startWeight; w <= endWeight; ++w)
				g_weight_to_secondIdx[w] = j;
		}

		/*compute linear upper bound (dual bound)*/
		DblMatrix ub_matr = { 0 };
		dm_new(&ub_matr, g_instance.n_items + 1, secondDim + 1);
		for (size_t i = 0; i <= g_instance.n_items; ++i) {
			for (size_t j = 0; j <= secondDim; ++j) {
				dm_set(&ub_matr, i, j, INFINITY);
			}
		}
		my_lin_relax(&g_instance, &ub_matr, g_instance.capacity + 1, weightRec, false);

		///*truncating a and b to get the dual bound*/ 
		//double dualBoundTime = 0;
		//CalculateDualbound(g_instance, weightRec, ub_matr, dualBoundTime);
		//cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;

		/*use the heuristic labeling algorithm to get better lower bound*/
		double heuLableTime = 0;
		unordered_map<int, ItemIndex> removedIndicesRec;
		vector<Stage> preStages;
		vector<double> sr3Duals;
		vector<unordered_set<int>> sr3s;
		LabelSettingHeuristic(ub_matr, heuLableTime, indicesRec, preStages, sr3Duals, sr3s, removedIndicesRec, finalSols);
		cout << "heuristic labeling get lower bound: " << g_bestLB << endl;
		cout << "heuristic labeling use time: " << heuLableTime << "s" << endl;

		//run the exact labeling algorithm
		multimap<double, KnapsackSol, greater<double>> preFinalSols;
		//initialize the bucket
		vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
		dominatedOldLabs.reserve(1e+6);
		int thisNonDominatedLabel = 0;
		if (!heuPricingFlag) {
			// truncating a and b to get the dual bound
			double dualBoundTime = 0;
			CalculateDualbound(g_instance, weightRec, ub_matr, sr3Duals, sr3s, preciseDuals, indicesRec, dualBoundTime);
			cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;

			auto exactDP_startTime = chrono::high_resolution_clock::now();
			vector<int> nonDominatedLabsRec(g_firstDim);

			MyLabel* initLab = new MyLabel();
			oldExtended->begin()->insert({ 0, initLab });
			++g_tolGeneratedLabel;

			//label extention and dominance
			int currItem = 0;				//record the current item to extend
			bool timeLimitFlag = false;
			for (int stage = 0; stage < g_firstDim; ++stage) {
				thisNonDominatedLabel = 0;
				for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
					if ((*oldExtended)[sndBktIdx].empty()) continue;
					auto ite = (*oldExtended)[sndBktIdx].begin();
					while (ite != (*oldExtended)[sndBktIdx].end()) {
						//check the if the old label can be discarded by CB
						bool keepOld = true;
						if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
							dominatedOldLabs.push_back(ite->second);
							keepOld = false;
						}
						else {
							++g_tolGeneratedLabel;
						}

						double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
							g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
						if (preWeight <= g_instance.capacity + EX) {//judge capacity
							//label extension
							MyLabel* tmpLab = new MyLabel(ite->second);
							LabelExtention(g_instance, ite->second, tmpLab, currItem);
							++g_tolGeneratedLabel;

							//completion bound to fathom label
							if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
								++g_CBFathomLabel;
								delete tmpLab; tmpLab = nullptr;
								//before insert the old label, do the dominance check
								MyLabel* oldLab = ite->second;
								ite = (*oldExtended)[sndBktIdx].erase(ite);
								if (keepOld) {
									bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
									if (dominanceFlag) {
										/*delete oldLab;
										oldLab = nullptr;*/
										dominatedOldLabs.push_back(oldLab);
										continue;
									}

									(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
								}
								continue;
							}

							//dominance check
							int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
							bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab);

							if (dominanceFlag) {
								delete tmpLab;
								++g_dominatedLabel;
							}
							else {
								//dominate other labels
								JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs);
								//add the non dominated label into the new extended bucket
								(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
							}
						}

						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
							}
							else {
								//save the current label to the new extended bucket
								(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
							}
						}
						////test maximum solution time
						//auto endTime = chrono::high_resolution_clock::now();
						//if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
						//	cout << "The maximum solution time of my DP is reached!" << endl;
						//	outPut << "The maximum solution time of my DP is reached!" << endl;
						//	timeLimitFlag = true;
						//	break;
						//}
					}
					//if (timeLimitFlag) break;
					thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
					//use the present best dual variable to update
					if (!(*newExtended)[sndBktIdx].empty()) {
						KnapsackSol preSol;
						preSol.bestLab = new MyLabel((*newExtended)[sndBktIdx].begin()->second);
						//get the item set
						if (preSol.bestLab->itemSet.empty())
							preSol.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
						auto tmpLab = preSol.bestLab;
						while (tmpLab->parentLab != nullptr) {
							preSol.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
							tmpLab = tmpLab->parentLab;
						}
						preFinalSols.insert({ (*newExtended)[sndBktIdx].begin()->first, preSol });
						g_bestLB = max(g_bestLB, (*newExtended)[sndBktIdx].begin()->first);
					}
				}
				auto endTime = chrono::high_resolution_clock::now();
				if (g_BPTime + chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
					cout << "The maximum solution time of my DP is reached!" << endl;
					outPut << "The maximum solution time of my DP is reached!" << endl;
					timeLimitFlag = true;
				}
				if (timeLimitFlag) {
					auto ite = (*oldExtended).begin();
					for (auto& e : (*newExtended)) {
						for (auto& t : e)
							ite->insert({ t.first, t.second });
						++ite;
					}
					break;
				}
				++currItem;
				delete oldExtended;
				oldExtended = newExtended;
				newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
				nonDominatedLabsRec[stage] = thisNonDominatedLabel;
			}
		}

		//record the best solution
		auto ite = --(*oldExtended).end();
		MyLabel* bestLab = nullptr;
		while (true) {
			//if (!ite->empty()) {
			//	bestLab = ite->begin()->second;
			//	ite->erase(ite->begin());
			//	break;
			//}
			//if (ite == (*oldExtended).begin())
			//	break;
			//--ite;
			if (!ite->empty()) {
				if (bestLab == nullptr)
					bestLab = ite->begin()->second;

				auto subIte = ite->begin();
				while (subIte != ite->end()) {
					/*if (1 - subIte->first >= 0)
						break;*/
					KnapsackSol preSol;
					preSol.bestLab = subIte->second;
					preFinalSols.insert({ subIte->second->tolProfit, preSol });
					subIte = ite->erase(subIte);
				}
			}
			if (ite == (*oldExtended).begin())
				break;
			--ite;
		}
		string bestItemSet = !finalSols.empty() ?
			finalSols.begin()->second.bestItemSet : "";
		if (bestLab != nullptr && bestLab->lastItem > 0 && bestLab->tolProfit > 1.0 + EX) {
			//get the item set
			vector<int> bestIS;
			MyLabel* tmpLab = bestLab;
			while (tmpLab->parentLab != nullptr) {
				bestIS.push_back(indicesRec[tmpLab->lastItem].index);
				tmpLab = tmpLab->parentLab;
			}
			sort(bestIS.begin(), bestIS.end());
			bestItemSet = JoinVector(bestIS);

			//store the best solution
			//KnapsackSol bestSol;
			//bestSol.bestItemSet = bestItemSet;
			//bestSol.bestLab = bestLab;
			//finalSols.insert({ bestLab->tolProfit, bestSol });
			g_bestLB = max(g_bestLB, preFinalSols.begin()->first);
			//get the item set for all the solutions
			for (auto& e : preFinalSols) {
				e.second.bestLab->itemSet.resize(g_instance.n_items, 0);
				tmpLab = e.second.bestLab;
				while (tmpLab->parentLab != nullptr) {
					e.second.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
					tmpLab = tmpLab->parentLab;
				}
			}
		}
		//merge sols
		for (auto& e : preFinalSols) {
			auto ite = finalSols.find(e.first);
			if (ite == finalSols.end())
				finalSols.insert(e);
			else
				delete e.second.bestLab;
		}
		//obtain valid columns
		auto tIte = finalSols.begin();
		if (!finalSols.empty())
			rcLB = min(0.0, 1.0 - tIte->first);
		while (tIte != finalSols.end()) {
			if (1.0 - tIte->first < -1e-6) {
				//if using upper bound of the dual variables, should test the real reduced cost
				if (g_upper_lower_dual == 0) {
					if (!TestPreciseReducedCost(tIte->second.bestLab, preciseDuals, sr3s, sr3Duals)) {
						delete tIte->second.bestLab;
						tIte = finalSols.erase(tIte);
						continue;
					}
				}
				//add new columns
				//newCols.emplace_back(tIte->second.bestLab);
				//tIte = finalSols.erase(tIte);
				//rcLB = min(rcLB, 1.0 - tIte->first);
				//duplication check
				if (find(allCols.begin(), allCols.end(), tIte->second.bestLab) == allCols.end()) {
					//add new columns
					newCols.emplace_back(tIte->second.bestLab);
					tIte = finalSols.erase(tIte);
				}
				else {
					delete tIte->second.bestLab;
					tIte = finalSols.erase(tIte);
					continue;
				}
				if (newCols.size() >= g_maxAddColNumPerCG)
					break;
				continue;
			}
			break;
		}

		//free space
		for (auto& e : finalSols)
			delete e.second.bestLab;
		for (auto& t : (*oldExtended)) {
			for (auto& e : t)
				delete e.second;
		}
		delete oldExtended;
		delete newExtended;
		g_dominatedLabel += dominatedOldLabs.size();
		for (auto& lab : dominatedOldLabs)
			delete lab;
		free(indicesRec);
		free(g_instance.a_ptr);
		free(g_instance.b_ptr);
		free(g_instance.p_ptr);
		free(g_instance.p_weight);
		g_instance.a_ptr = g_instance.b_ptr = g_instance.p_ptr = g_instance.p_weight = nullptr;
		dm_free(&ub_matr);

		cout << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
		cout << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
		cout << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
		cout << "The number of dominated labels is: " << g_dominatedLabel << endl;
		g_tolGeneratedLabel = 0;
		g_nonDominatedLabel = 0;
		g_CBFathomLabel = 0;
		g_dominatedLabel = 0;
		auto endTime = chrono::high_resolution_clock::now();
		preTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
		return !newCols.empty();
	}
	catch (const exception& e) {
		g_errorControl = true;
		cerr << "Error in LabelSettingSolveSubproblem: " << e.what() << endl;
		outPut << "Error in LabelSettingSolveSubproblem: " << e.what() << endl;
		return !newCols.empty();
	}
}

//adjust instances according to the preStages
void removeAndShift(
	double* a,
	double* b,
	double* p,
	double* p_w,
	ItemIndex*& indiciesRec,
	unordered_map<int, ItemIndex>& removedIndicesRec,
	int n,
	int indexToRemove
) {
	if (indexToRemove < 0 || indexToRemove >= n) {
		std::cerr << "Invalid indexToRemove index!" << std::endl;
		exit(-1);
	}

	double temp_a = a[indexToRemove];  // save the value of to be removed item
	double temp_b = b[indexToRemove];
	double temp_p = p[indexToRemove];
	double temp_pw = p_w[indexToRemove];
	ItemIndex temp_indice = indiciesRec[indexToRemove];

	// remaining item move forward
	for (int i = indexToRemove; i < n - 1; ++i) {
		a[i] = a[i + 1];
		b[i] = b[i + 1];
		p[i] = p[i + 1];
		p_w[i] = p_w[i + 1];
		indiciesRec[i] = indiciesRec[i + 1];
	}

	// the deleted item removed to the tail of the instances
	a[n - 1] = temp_a;
	b[n - 1] = temp_b;
	p[n - 1] = temp_p;
	p_w[n - 1] = temp_pw;
	indiciesRec[n - 1] = temp_indice;
	temp_indice.pos = n - 1;
	removedIndicesRec.insert({ temp_indice.index, temp_indice });
}
void AdjustInstancesPerStages(vector<Stage>& preStages, ItemIndex*& indicesRec, unordered_map<int, ItemIndex>& removedIndicesRec) {
	//remove all the items in preStages
	unordered_map<int, int> orgIdx_sortIdx;
	for (int i = 0; i < g_instance.n_items; ++i)
		orgIdx_sortIdx[indicesRec[i].index] = i;
	vector<int> itemsIdxToRmv;
	for (auto& sg : preStages) {
		for (auto& e : sg.stageCandidates)
			itemsIdxToRmv.push_back(orgIdx_sortIdx[e.first]);
	}
	sort(itemsIdxToRmv.begin(), itemsIdxToRmv.end());//first remove the item with larger index
	for (int i = itemsIdxToRmv.size() - 1; i >= 0; --i) {
		//remove one item
		removeAndShift(g_instance.a_ptr, g_instance.b_ptr, g_instance.p_ptr, g_instance.p_weight, indicesRec, removedIndicesRec, g_instance.n_items, itemsIdxToRmv[i]);
		--g_instance.n_items;
	}
}

//use the label algorithm as the subproblem of the bin packing problem, used for child nodes
bool LabelSettingSolveSubproblem_Child(
	vector<Bin*>& newCols,
	ofstream& outPut,
	const vector<Bin*>& allCols,
	vector<Stage>& preStages,
	double& preTime,
	vector<double>& preciseDuals,
	vector<double>& SR3Duals,
	vector<unordered_set<int>>& SR3s,
	double& rcLB,
	unordered_map<int, BranchPattern>& branchPatterns,
	bool heuPricingFlag
) {
	try {
		auto startTime = chrono::high_resolution_clock::now();
		//sort the instance according to p/(a+sqrt(b))
		ItemIndex* indicesRec = nullptr;
		sort_instance_by_p_ab(g_instance, indicesRec);

		/*adjust instances according to the preStages*/
		unordered_map<int, ItemIndex> removedIndicesRec;
		AdjustInstancesPerStages(preStages, indicesRec, removedIndicesRec);


		//specify the dimension of the bucket
		g_firstDim = g_instance.n_items / 1;		//the first dimension of the bucket
		secondDim = g_instance.capacity / g_secondDimDiv;       //the second dimension of the bucket

		/*use heuristic to get primal bound*/
		//get the initial lower bound
		multimap<double, KnapsackSol, greater<double>> finalSols;
		g_bestLB = 0;
		GetInitialLowerBound(g_bestLB);
		//TS heuristic
		vector<int> ts_sol_items(g_tolItemNum);
		double ts_sol = heuristic_solution(&g_instance, &ts_sol_items[0]);
		double tsSol_sum_a = 0;
		double tsSol_sum_b = 0;
		vector<int> map_ts_sol_items(g_tolItemNum);
		for (int i = 0; i < g_tolItemNum; ++i) {
			if (ts_sol_items[i]) {
				tsSol_sum_a += g_instance.a_ptr[i];
				tsSol_sum_b += g_instance.b_ptr[i];
				map_ts_sol_items[indicesRec[i].index] = 1;
			}
		}
		MyLabel* tsSolLab = new MyLabel();
		tsSolLab->tolProfit = ts_sol;
		tsSolLab->itemSet = map_ts_sol_items;
		tsSolLab->tolWeight = tsSol_sum_a + g_instance.rho * sqrt(tsSol_sum_b);
		KnapsackSol preSol;
		preSol.bestLab = tsSolLab;
		finalSols.insert({ ts_sol, preSol });
		//the best primal bound is at least equal to 1
		if (preStages.empty())
			g_bestLB = max(max(g_bestLB, 1.0), ts_sol);

		/*construct bucket index*/
		int itemGap = g_instance.n_items / (g_firstDim > 0 ? g_firstDim : 1);	//the item gap for each item bucket
		g_item_to_firstIdx.resize(g_instance.n_items + 1, 0);
		for (int j = 0; j <= g_firstDim; ++j) {
			int startItem = j * itemGap;
			int endItem = (j == g_firstDim - 1) ? g_instance.n_items : (j + 1) * itemGap - 1;
			endItem = min((size_t)endItem, g_instance.n_items);
			for (int w = startItem; w <= endItem; ++w)
				g_item_to_firstIdx[w] = j;
		}
		vector<double> weightRec(secondDim + 1, 0);
		int weightGap = g_instance.capacity / secondDim;	//the capacity gap for each capacity bucket
		g_weight_to_secondIdx.resize(g_instance.capacity + 1, 0);
		for (int j = 0; j <= secondDim; ++j) {
			int startWeight = j * weightGap;
			int endWeight = (j == secondDim - 1) ? g_instance.capacity : (j + 1) * weightGap - 1;
			endWeight = min(endWeight * 1.0, g_instance.capacity);
			weightRec[j] = endWeight;
			for (int w = startWeight; w <= endWeight; ++w)
				g_weight_to_secondIdx[w] = j;
		}

		/*compute linear upper bound (dual bound)*/
		DblMatrix ub_matr = { 0 };
		if (g_instance.n_items > 0) {
			dm_new(&ub_matr, g_instance.n_items + 1, secondDim + 1);
			for (size_t i = 0; i <= g_instance.n_items; ++i) {
				for (size_t j = 0; j <= secondDim; ++j) {
					dm_set(&ub_matr, i, j, INFINITY);
				}
			}
			my_lin_relax(&g_instance, &ub_matr, g_instance.capacity + 1, weightRec, false);
		}

		///*truncating a and b to get the dual bound*/ 
		//double dualBoundTime = 0;
		//CalculateDualbound(g_instance, weightRec, ub_matr, dualBoundTime, binNumDual);
		//cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;

		/*use the heuristic labeling algorithm to get better lower bound*/
		double heuLableTime = 0;
		LabelSettingHeuristic(ub_matr, heuLableTime, indicesRec, preStages, SR3Duals, SR3s, removedIndicesRec, finalSols);
		cout << "heuristic labeling get lower bound: " << g_bestLB << endl;
		cout << "heuristic labeling use time: " << heuLableTime << "s" << endl;

		//run the exact labeling algorithm
		multimap<double, KnapsackSol, greater<double>> preFinalSols;
		//initialize the bucket
		vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels
		dominatedOldLabs.reserve(1e+6);
		int thisNonDominatedLabel = 0;
		if (!heuPricingFlag) {
			// truncating a and b to get the dual bound
			if (g_instance.n_items > 0) {
				double dualBoundTime = 0;
				CalculateDualbound(g_instance, weightRec, ub_matr, SR3Duals, SR3s, preciseDuals, indicesRec, dualBoundTime);
				cout << "The time of calculating dual bound is: " << dualBoundTime << "s" << endl;
			}
			auto exactDP_startTime = chrono::high_resolution_clock::now();
			vector<int> nonDominatedLabsRec(g_firstDim);

			MyLabel* initLab = new MyLabel();
			oldExtended->begin()->insert({ 0, initLab });
			++g_tolGeneratedLabel;

			//label extention and dominance
			int currItem = 0;				//record the current item to extend
			bool timeLimitFlag = false;
			for (int stage = 0; stage < g_firstDim; ++stage) {
				thisNonDominatedLabel = 0;
				for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
					if ((*oldExtended)[sndBktIdx].empty()) continue;
					auto ite = (*oldExtended)[sndBktIdx].begin();
					while (ite != (*oldExtended)[sndBktIdx].end()) {
						//check the if the old label can be discarded by CB
						bool keepOld = true;
						if (CompletionBound(ite->second, currItem + 1, &ub_matr)) {
							dominatedOldLabs.push_back(ite->second);
							keepOld = false;
						}
						else {
							++g_tolGeneratedLabel;
						}

						double preWeight = ite->second->sum_a + g_instance.a_ptr[currItem] +
							g_instance.rho * sqrt(ite->second->sum_b + g_instance.b_ptr[currItem]);
						if (preWeight <= g_instance.capacity + EX) {//judge capacity
							//label extension
							MyLabel* tmpLab = new MyLabel(ite->second);
							LabelExtention(g_instance, ite->second, tmpLab, currItem);
							++g_tolGeneratedLabel;

							//completion bound to fathom label
							if (CompletionBound(tmpLab, currItem + 1, &ub_matr)) {
								++g_CBFathomLabel;
								delete tmpLab; tmpLab = nullptr;
								//before insert the old label, do the dominance check
								MyLabel* oldLab = ite->second;
								ite = (*oldExtended)[sndBktIdx].erase(ite);
								if (keepOld) {
									bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
									if (dominanceFlag) {
										/*delete oldLab;
										oldLab = nullptr;*/
										dominatedOldLabs.push_back(oldLab);
										continue;
									}

									(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
								}
								continue;
							}

							//dominance check
							int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
							bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab);

							if (dominanceFlag) {
								delete tmpLab;
								++g_dominatedLabel;
							}
							else {
								//dominate other labels
								JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs);
								//add the non dominated label into the new extended bucket
								(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
							}
						}

						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false);
							if (dominanceFlag) {
								/*delete oldLab;
								oldLab = nullptr;*/
								dominatedOldLabs.push_back(oldLab);
							}
							else {
								//save the current label to the new extended bucket
								(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
							}
						}
						////test maximum solution time
						//auto endTime = chrono::high_resolution_clock::now();
						//if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
						//	cout << "The maximum solution time of my DP is reached!" << endl;
						//	outPut << "The maximum solution time of my DP is reached!" << endl;
						//	timeLimitFlag = true;
						//	break;
						//}
					}
					//if (timeLimitFlag) break;
					thisNonDominatedLabel += (*newExtended)[sndBktIdx].size();
					//use the present best dual variable to update
					if (!(*newExtended)[sndBktIdx].empty()) {
						KnapsackSol preSol;
						preSol.bestLab = new MyLabel((*newExtended)[sndBktIdx].begin()->second);
						//get the item set
						if (preSol.bestLab->itemSet.empty())
							preSol.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
						auto tmpLab = preSol.bestLab;
						while (tmpLab->parentLab != nullptr) {
							preSol.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
							tmpLab = tmpLab->parentLab;
						}
						preFinalSols.insert({ (*newExtended)[sndBktIdx].begin()->first, preSol });
						g_bestLB = max(g_bestLB, (*newExtended)[sndBktIdx].begin()->first);
					}
				}
				auto endTime = chrono::high_resolution_clock::now();
				if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 >= MAXMUMSOLTIME) {
					cout << "The maximum solution time of my DP is reached!" << endl;
					outPut << "The maximum solution time of my DP is reached!" << endl;
					timeLimitFlag = true;
				}
				if (timeLimitFlag) {
					auto ite = (*oldExtended).begin();
					for (auto& e : (*newExtended)) {
						for (auto& t : e)
							ite->insert({ t.first, t.second });
						++ite;
					}
					break;
				}
				++currItem;
				delete oldExtended;
				oldExtended = newExtended;
				newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
				nonDominatedLabsRec[stage] = thisNonDominatedLabel;
			}

			//if Stages is not empty, do the final extensions
			if (!preStages.empty()) {//add new valid labels into oldExtended
				ExtendNDlabsToStages(oldExtended, dominatedOldLabs, removedIndicesRec, preStages, SR3Duals);
			}
		}

		//record the best solution
		auto ite = --(*oldExtended).end();
		MyLabel* bestLab = nullptr;
		while (true) {
			//if (!ite->empty()) {
			//	bestLab = ite->begin()->second;
			//	ite->erase(ite->begin());
			//	break;
			//}
			//if (ite == (*oldExtended).begin())
			//	break;
			//--ite;
			if (!ite->empty()) {
				if (bestLab == nullptr)
					bestLab = ite->begin()->second;

				auto subIte = ite->begin();
				while (subIte != ite->end()) {
					/*if (1 - subIte->first >= 0)
						break;*/
					KnapsackSol preSol;
					preSol.bestLab = subIte->second;
					preFinalSols.insert({ subIte->second->tolProfit, preSol });
					subIte = ite->erase(subIte);
				}
			}
			if (ite == (*oldExtended).begin())
				break;
			--ite;
		}
		string bestItemSet = !finalSols.empty() ?
			finalSols.begin()->second.bestItemSet : "";
		if (bestLab != nullptr && bestLab->lastItem > 0 && bestLab->tolProfit > 1.0 + EX) {
			//get the item set
			vector<int> bestIS;
			MyLabel* tmpLab = bestLab;
			while (tmpLab->parentLab != nullptr) {
				bestIS.push_back(indicesRec[tmpLab->lastItem].index);
				tmpLab = tmpLab->parentLab;
			}
			sort(bestIS.begin(), bestIS.end());
			bestItemSet = JoinVector(bestIS);

			//store the best solution
			//KnapsackSol bestSol;
			//bestSol.bestItemSet = bestItemSet;
			//bestSol.bestLab = bestLab;
			//finalSols.insert({ bestLab->tolProfit, bestSol });
			g_bestLB = max(g_bestLB, preFinalSols.begin()->first);
			//get the item set for all the solutions
			for (auto& e : preFinalSols) {
				if (e.second.bestLab->itemSet.empty())
					e.second.bestLab->itemSet.resize(g_instance.n_items + removedIndicesRec.size());
				tmpLab = e.second.bestLab;
				while (tmpLab->parentLab != nullptr) {
					e.second.bestLab->itemSet[indicesRec[tmpLab->lastItem].index] = 1;
					tmpLab = tmpLab->parentLab;
				}
			}
		}
		//merge sols
		for (auto& e : preFinalSols) {
			auto ite = finalSols.find(e.first);
			if (ite == finalSols.end())
				finalSols.insert(e);
			else
				delete e.second.bestLab;
		}
		//obtain valid columns
		auto tIte = finalSols.begin();
		if (!finalSols.empty())
			rcLB = min(0.0, 1.0 - tIte->first);
		while (tIte != finalSols.end()) {
			if (1.0 - tIte->first < -1e-6) {
				//if using lower bound of the dual variables, should test the real reduced cost
				if (g_upper_lower_dual == 0) {
					if (!TestPreciseReducedCost(tIte->second.bestLab, preciseDuals, SR3s, SR3Duals)) {
						delete tIte->second.bestLab;
						tIte = finalSols.erase(tIte);
						continue;
					}
				}
				//check if the column violate the branch patterns
				bool vioFlag = false;
				for (auto& pt : branchPatterns) {
					if (
						tIte->second.bestLab->itemSet[pt.first] &&
						(tIte->second.bestLab->tolWeight < pt.second.branchRange.first - EX ||
							tIte->second.bestLab->tolWeight > pt.second.branchRange.second + EX)
						) {
						vioFlag = true;
						break; //skip the corrent column
					}
				}
				if (vioFlag) {
					delete tIte->second.bestLab;
					tIte = finalSols.erase(tIte);
					continue;
				}
				//add new columns
				/*newCols.emplace_back(tIte->second.bestLab);
				tIte = finalSols.erase(tIte);*/
				//rcLB = min(rcLB, 1.0 - tIte->first);
				//duplication check
				if (find(allCols.begin(), allCols.end(), tIte->second.bestLab) == allCols.end()) {
					//add new columns
					newCols.emplace_back(tIte->second.bestLab);
					tIte = finalSols.erase(tIte);
				}
				else {
					delete tIte->second.bestLab;
					tIte = finalSols.erase(tIte);
					continue;
				}

				if (newCols.size() >= g_maxAddColNumPerCG)
					break;
				continue;
			}
			break;
		}

		//free space
		for (auto& e : finalSols)
			delete e.second.bestLab;
		for (auto& t : (*oldExtended)) {
			for (auto& e : t)
				delete e.second;
		}
		delete oldExtended;
		delete newExtended;
		g_dominatedLabel += dominatedOldLabs.size();
		for (auto& lab : dominatedOldLabs)
			delete lab;
		free(indicesRec);
		free(g_instance.a_ptr);
		free(g_instance.b_ptr);
		free(g_instance.p_ptr);
		free(g_instance.p_weight);
		g_instance.a_ptr = g_instance.b_ptr = g_instance.p_ptr = g_instance.p_weight = nullptr;
		dm_free(&ub_matr);

		cout << "The total number of generated labels is: " << g_tolGeneratedLabel << endl;
		cout << "The number of non-dominated labels is: " << g_nonDominatedLabel + thisNonDominatedLabel << endl;
		cout << "The number of CB fathomed labels is: " << g_CBFathomLabel << endl;
		cout << "The number of dominated labels is: " << g_dominatedLabel << endl;
		g_tolGeneratedLabel = 0;
		g_nonDominatedLabel = 0;
		g_CBFathomLabel = 0;
		g_dominatedLabel = 0;
		auto endTime = chrono::high_resolution_clock::now();
		preTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
		return !newCols.empty();
	}
	catch (const exception& e) {
		g_errorControl = true;
		cerr << "Error in LabelSettingSolveSubproblem: " << e.what() << endl;
		outPut << "Error in LabelSettingSolveSubproblem: " << e.what() << endl;
		return !newCols.empty();
	}
}