#include "labelSetting.h"
//the maximum time that can spent on calcalation the dual bound

using namespace std;

int g_dualBoundRoundNum = 1;			//the roundNum when executing the dual bound calculation
int g_maximumTimeDualBound = 300;		//the time limitation for executing the dual bound calculation

// ascending sort
int ReverseSort(const void* a, const void* b) {
	double diff = ((ItemIndex*)a)->sum_p_ab - ((ItemIndex*)b)->sum_p_ab;
	return (diff > 0) - (diff < 0); // ·µ»Ø1¡¢0»ò-1
}

// modified_weight, descending sort
int compare_items_by_modified_weight(const void* a, const void* b) {
	const Item* item1 = (const Item*)a;
	const Item* item2 = (const Item*)b;

	double w1 = item1->a_val + g_instance.rho * sqrt(g_smallest_bDiva * item1->a_val);
	double w2 = item2->a_val + g_instance.rho * sqrt(g_smallest_bDiva * item2->a_val);

	w1 = item1->p_val / w1;
	w2 = item2->p_val / w2;

	if (w1 < w2 - EX) return 1;
	if (w1 > w2 + EX) return -1;
	return 0;
}

//build the projection from the one weight to the index of the second bucket
void build_weight_to_secondIdx(int* weight_to_secondIdx, int capacity, int secondDim) {
	// initialize the array to 0
	memset(weight_to_secondIdx, 0, (capacity + 1) * sizeof(int));
	int weightGap = capacity / (secondDim - 1);

	for (int j = 0; j <= secondDim; ++j) {
		int startWeight = j * weightGap;
		int endWeight;
		if (j == secondDim - 1) {
			endWeight = capacity;
		}
		else {
			endWeight = (j + 1) * weightGap - 1;
		}
		if (endWeight > capacity) {
			endWeight = capacity;
		}
		for (int w = startWeight; w <= endWeight; ++w) {
			weight_to_secondIdx[w] = j;
		}
	}
}
int my_lin_relax(const Instance* RESTRICT const inst,
	DblMatrix* RESTRICT const output,
	int binCapacity,
	const vector<double>& weightRec,
	bool dualFlag) {

	////build the projection from the one weight to the index of the second bucket
	//int* weight_to_secondIdx = (int*)malloc(binCapacity * sizeof(int));
	//build_weight_to_secondIdx(weight_to_secondIdx, binCapacity - 1, output->rows);

	ItemList item_list = { 0 };
	run_or_fail(item_list_alloc(&item_list, inst, FullWeight));
	if (!dualFlag)
		invert_item_list(&item_list);

	size_t n_rows = output->cols;
	item_list.n_items = 0;
	//size_t n_cols = output->rows;
	//size_t n_cols = binCapacity;
	size_t n_cols = weightRec.size() - 1;

	for (size_t i = 1; i < n_rows; ++i) {
		++item_list.n_items;
		item_list_sort(&item_list, AWeight);
		////sort according to modified weight
		//qsort(item_list.item_ptr, item_list.n_items, sizeof(Item), compare_items_by_modified_weight);

		double last_w = 0;
		double last_p = 0;
		size_t start_item = 0;
		bool complete = false;
		for (size_t j = 0; j < n_cols; ++j) {
			//double inW = j + 1;
			double w = weightRec[j + 1];
			double curr_w = last_w;
			//double curr_w = last_w + inst->rho * sqrt(g_smallest_bDiva * last_w);
			double curr_p = last_p;
			bool run = true;
			if (!complete) {
				for (size_t k = start_item; k < i && run; ++k) {
					Item item = item_list.item_ptr[k];
					double avail = w - curr_w;
					double inW = item.a_val;
					//double inW = item.a_val + inst->rho * sqrt(g_smallest_bDiva * item.a_val);
					double p = item.p_val;
					//if (avail < item.a_val) {
					if (avail < inW) {
						//double frac = avail / item.a_val;
						double frac = avail / inW;
						inW *= frac;
						p *= frac;
						run = false;
						start_item = k;
					}
					else {
						last_w += inW;
						last_p += p;
					}
					curr_w += inW;
					curr_p += p;
				}
				if (run)
					complete = true;
			}
			////map j to the index of the second bucket
			//int secondIdx = weight_to_secondIdx[j];
			//if (j == n_cols - 1 || weight_to_secondIdx[j + 1] > secondIdx) {
			//	dm_set(output, i, secondIdx, curr_p);
			//}
			dm_set(output, i, j, curr_p);
		}
	}

	item_list_free(&item_list);
	//free(weight_to_secondIdx);
	return SUCCESS;
}

//floor to k decimal places (return new value)
double roundToKDecimal(double value, int k) {
	double factor = std::pow(10.0, k);
	return std::floor(value * factor) / factor;
}

//convert instance by truncating a and b and retain some decimal
void ConvertInstance(
	Instance& g_Instance,
	Instance& thisInstance,
	int roundNum,
	vector<double>& preciseDuals,
	const ItemIndex* indicesRec
) {
	thisInstance.n_items = g_Instance.n_items;
	thisInstance.capacity = g_Instance.capacity;
	thisInstance.rho = g_Instance.rho;
	thisInstance.a_ptr = new double[thisInstance.n_items];
	thisInstance.b_ptr = new double[thisInstance.n_items];
	thisInstance.p_ptr = new double[thisInstance.n_items];
	thisInstance.p_weight = new double[thisInstance.n_items];

	double factor = 1;
	if (g_dual_bound_parameter == 1)//use variable parameters
		factor = std::pow(10.0, g_dualBoundRoundNum + 1);
	for (size_t i = 0; i < g_Instance.n_items; ++i) {
		thisInstance.a_ptr[i] = roundToKDecimal(g_Instance.a_ptr[i], roundNum);
		thisInstance.b_ptr[i] = roundToKDecimal(g_Instance.b_ptr[i], roundNum);
		if (g_dual_bound_parameter == 0)//use fixed parameters
			thisInstance.p_ptr[i] = g_Instance.p_ptr[i];
		else {//use variable parameters
			if (factor < g_controlNum) {
				//use upper bound of the dual variables
				thisInstance.p_ptr[i] = ceil(factor * preciseDuals[indicesRec[i].index]) / factor;
			}
			else
				thisInstance.p_ptr[i] = g_Instance.p_ptr[i];
		}

		thisInstance.p_weight[i] = g_Instance.p_weight[i];
	}
}

//truncate a and b to speed up the algorithm and get dual bound 
void CalculateDualbound(
	Instance& g_Instance,
	const vector<double>& weightRec,
	DblMatrix& ub_matr,
	vector<double> SR3Duals,
	vector<unordered_set<int>> SR3s,
	vector<double>& preciseDuals,
	const ItemIndex* indicesRec,
	double& dualBoundTime
) {
	auto startTime = chrono::high_resolution_clock::now();

	//convert instance by truncating a and b and retain some decimal
	Instance thisInstance;
	if (g_dual_bound_parameter == 0) {//use fixed parameters
		int roundNum = 1;	//the number of decimal to keep
		ConvertInstance(g_Instance, thisInstance, roundNum, preciseDuals, indicesRec);
	}
	else {//use variable parameters
		ConvertInstance(g_Instance, thisInstance, g_dualBoundRoundNum, preciseDuals, indicesRec);
	}

	//reverse the order of the instance
	std::reverse(thisInstance.a_ptr, thisInstance.a_ptr + thisInstance.n_items);
	std::reverse(thisInstance.b_ptr, thisInstance.b_ptr + thisInstance.n_items);
	std::reverse(thisInstance.p_ptr, thisInstance.p_ptr + thisInstance.n_items);
	std::reverse(thisInstance.p_weight, thisInstance.p_weight + thisInstance.n_items);

	//compute linear upper bound (dual bound)
	DblMatrix pre_ub_matr = { 0 };
	dm_new(&pre_ub_matr, thisInstance.n_items + 1, secondDim + 1);
	for (size_t i = 0; i < thisInstance.n_items; ++i)
		for (size_t j = 0; j < secondDim; ++j)
			dm_set(&pre_ub_matr, i, j, INFINITY);
	my_lin_relax(&thisInstance, &pre_ub_matr, thisInstance.capacity + 1, weightRec, true);

	//initialize the bucket
	vector<multimap<double, MyLabel*, greater<double>>>* oldExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<multimap<double, MyLabel*, greater<double>>>* newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
	vector<MyLabel*> dominatedOldLabs;	//record the dominated old labels

	MyLabel* initLab = new MyLabel();
	oldExtended->begin()->insert({ 0, initLab });

	int improveNum = 0;
	//label extention and dominance
	int currItem = 0;				//record the current item to extend
	for (int stage = 0; stage < g_firstDim; ++stage) {
		for (int sndBktIdx = 0; sndBktIdx <= secondDim; ++sndBktIdx) {
			if ((*oldExtended)[sndBktIdx].empty()) continue;
			int i = stage + 1;
			int j = sndBktIdx - 1; //remain capacity
			////use the best labels in the oldExtended to update ub_matr
			//if (!(*oldExtended)[sndBktIdx].empty() &&
			//	j >= 0 &&
			//	(*oldExtended)[sndBktIdx].begin()->first < dm_get(&ub_matr, i, j) - EX) {
			//	/*cout << "Old lable hit dual bound update " << i << "-" << j << "  " <<
			//		dm_get(&ub_matr, i, j) << " V.S. " <<
			//		(*oldExtended)[sndBktIdx].begin()->first << endl;*/
			//	dm_set(&ub_matr, i, j, (*oldExtended)[sndBktIdx].begin()->first);
			//	++improveNum;
			//}

			auto ite = (*oldExtended)[sndBktIdx].begin();
			while (ite != (*oldExtended)[sndBktIdx].end()) {
				//check the if the old label can be discarded by CB
				bool keepOld = true;
				if (CompletionBound(ite->second, currItem + 1, &pre_ub_matr)) {
					//dominatedOldLabs.push_back(ite->second);
					keepOld = false;
				}

				double preWeight = ite->second->sum_a + thisInstance.a_ptr[currItem] +
					thisInstance.rho * sqrt(ite->second->sum_b + thisInstance.b_ptr[currItem]);
				if (preWeight <= thisInstance.capacity + EX) {//judge capacity
					//label extension
					MyLabel* tmpLab = new MyLabel(ite->second);
					LabelExtention(thisInstance, ite->second, tmpLab, currItem);
					++g_tolGeneratedLabel;

					////use each generated label to update completion bound
					//int preJ = g_weight_to_secondIdx[(size_t)ceil(tmpLab->tolWeight)] - 1;
					//if (preJ >= 0 &&
					//	tmpLab->tolProfit < dm_get(&ub_matr, i, preJ) - EX) {
					//	/*cout << "New generated label hit dual bound update " << i << "-" << preJ << "  " <<
					//		dm_get(&ub_matr, i, preJ) << " V.S. " <<
					//		tmpLab->tolProfit << endl;*/
					//	dm_set(&ub_matr, i, preJ, tmpLab->tolProfit);
					//	++improveNum;
					//}

					//completion bound to fathom label
					if (CompletionBound(tmpLab, currItem + 1, &pre_ub_matr)) {
						delete tmpLab; tmpLab = nullptr;
						//before insert the old label, do the dominance check
						MyLabel* oldLab = ite->second;
						ite = (*oldExtended)[sndBktIdx].erase(ite);
						if (keepOld) {
							bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, false);
							if (dominanceFlag) {
								delete oldLab;
								oldLab = nullptr;
								//dominatedOldLabs.push_back(oldLab);
								continue;
							}
							(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
						}
						else {
							delete oldLab;
							oldLab = nullptr;
						}
						continue;
					}

					//dominance check
					int preSecondIdx = g_weight_to_secondIdx[(size_t)ceil(preWeight)];
					bool dominanceFlag = DominanceLogic(preSecondIdx, newExtended, oldExtended, tmpLab, true, false);

					if (dominanceFlag) {
						delete tmpLab;
						tmpLab = nullptr;
					}
					else {
						//dominate other labels
						JgeDominance((*newExtended), (*oldExtended), tmpLab, dominatedOldLabs, true, false, true);
						//add the non dominated label into the new extended bucket
						(*newExtended)[preSecondIdx].insert({ tmpLab->tolProfit, tmpLab });
					}
				}

				//before insert the old label, do the dominance check
				MyLabel* oldLab = ite->second;
				ite = (*oldExtended)[sndBktIdx].erase(ite);
				if (keepOld) {
					bool dominanceFlag = DominanceLogic(sndBktIdx, newExtended, oldExtended, oldLab, false, false);
					if (dominanceFlag) {
						delete oldLab;
						oldLab = nullptr;
						//dominatedOldLabs.push_back(oldLab);
						continue;
					}
					else {
						//save the current label to the new extended bucket
						(*newExtended)[sndBktIdx].insert({ oldLab->tolProfit, oldLab });
					}
				}
				else {
					delete oldLab;
					oldLab = nullptr;
				}
			}
			//use the best labels in the newExtended to update ub_matr
			if (!(*newExtended)[sndBktIdx].empty() &&
				j >= 0 &&
				(*newExtended)[sndBktIdx].begin()->first < dm_get(&ub_matr, i, j) - EX) {
				/*cout << "Hit dual bound update " << i << "-" << j << "  " <<
					dm_get(&ub_matr, i, j) << " V.S. " <<
					(*newExtended)[sndBktIdx].begin()->first << endl;*/
				dm_set(&ub_matr, i, j, (*newExtended)[sndBktIdx].begin()->first);
				++improveNum;
			}
			//dualbound calculation time limitation
			auto endTime = chrono::high_resolution_clock::now();
			if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 > g_maximumTimeDualBound)
				break;
		}

		delete oldExtended;
		oldExtended = newExtended;
		newExtended = new vector<multimap<double, MyLabel*, greater<double>>>(secondDim + 1, multimap<double, MyLabel*, greater<double>>());
		++currItem;
		//dualbound calculation time limitation
		auto endTime = chrono::high_resolution_clock::now();
		if (chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 > g_maximumTimeDualBound)
			break;
	}
	cout << "Total improve number is " << improveNum << endl;

	//free space
	for (auto& t : (*oldExtended)) {
		for (auto& e : t) {
			delete e.second;
			e.second = nullptr;
		}
	}
	delete oldExtended;
	delete newExtended;
	for (auto& t : dominatedOldLabs) {
		delete t;
		t = nullptr;
	}
	//free thisInstance
	delete[] thisInstance.a_ptr;
	delete[] thisInstance.b_ptr;
	delete[] thisInstance.p_ptr;
	delete[] thisInstance.p_weight;
	thisInstance.a_ptr = thisInstance.b_ptr = thisInstance.p_ptr = thisInstance.p_weight = nullptr;
	dm_free(&pre_ub_matr);
	auto endTime = chrono::high_resolution_clock::now();
	dualBoundTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
}