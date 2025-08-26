#include "labelSetting.h"

using namespace std;

int g_populationSize = 100;
int g_populationMaxSize = 110;

//initial solution
void GenerateInitialSol(vector<Bin*>& initSol, Instance& preInstances) {
    Bin* preBin = new Bin();
    preBin->itemSet.resize(preInstances.n_items);
    for (int i = 0; i < preInstances.n_items; ++i) {
        //judge the capacity constraints
        double tmpWeight = preBin->sum_a + preInstances.a_ptr[i] +
            preInstances.rho * sqrt(preBin->sum_b + preInstances.b_ptr[i]);

        if (tmpWeight <= preInstances.capacity + EX) {//add item to the current bin
            preBin->itemSet[i] = 1;
            preBin->sum_a += preInstances.a_ptr[i];
            preBin->sum_b += preInstances.b_ptr[i];
            preBin->tolWeight = tmpWeight;
        }
        else {//open a new bin
            initSol.emplace_back(preBin);
            preBin = new Bin();
            preBin->itemSet.resize(preInstances.n_items);
            --i;
        }
    }
    if (preBin->sum_a > 0)
        initSol.emplace_back(preBin);
}

//test bin capacity
bool TestBinCapacity(vector<int>& itemSet, Instance& preInstances) {
    double preWeight = 0;
    double sum_a = 0;
    double sum_b = 0;
    for (int i = 0; i < itemSet.size(); ++i) {
        if (itemSet[i]) {
            sum_a += preInstances.a_ptr[i];
            sum_b += preInstances.b_ptr[i];
            preWeight = sum_a + preInstances.rho * sqrt(sum_b);
            if (preWeight > preInstances.capacity) {
                cerr << "Bin capacity violated!" << endl;
                return false;
            }
        }
    }
    return true;
}

//LS
void LocalSearch(vector<Bin*>& sol, Instance& inst) {
    int i = rand() % sol.size();
    int j = rand() % (sol.size() - 1);
    if (j >= i)
        ++j;
    Bin* bin1 = sol[i];
    Bin* bin2 = sol[j];

    //swap
    int randStartIdx = rand() % inst.n_items;
    for (int k = randStartIdx; k < inst.n_items; ++k) {
        if (bin1->itemSet[k] && !bin2->itemSet[k]) {
            double newWeight2 = bin2->sum_a + inst.a_ptr[k] + inst.rho * sqrt(bin2->sum_b + inst.b_ptr[k]);
            if (newWeight2 <= inst.capacity) {
                bin1->itemSet[k] = 0;
                bin2->itemSet[k] = 1;
                bin1->sum_a -= inst.a_ptr[k];
                bin1->sum_b -= inst.b_ptr[k];
                bin1->tolWeight = bin1->sum_a + inst.rho * sqrt(bin1->sum_b);
                bin2->sum_a += inst.a_ptr[k];
                bin2->sum_b += inst.b_ptr[k];
                bin2->tolWeight = newWeight2;
            }
        }
    }

    //2-opt
    int n = sol.size();
    bool improved = true;
    while (improved) {
        improved = false;
        int i = rand() % sol.size();
        int j = rand() % (sol.size() - 1);
        if (j >= i)
            ++j;
        Bin* bin1 = sol[i];
        Bin* bin2 = sol[j];
        // 获取 bin 中 item 索引集合
        vector<int> items1, items2;
        for (int k = 0; k < inst.n_items; ++k) {
            if (bin1->itemSet[k]) items1.push_back(k);
            if (bin2->itemSet[k]) items2.push_back(k);
        }
        // 遍历所有可能交换的 pairre
        for (int idx1 : items1) {
            int cnt = 0;
            for (int idx2 : items2) {
                // 尝试交换 idx1 和 idx2
                double new_a1 = bin1->sum_a - inst.a_ptr[idx1] + inst.a_ptr[idx2];
                double new_b1 = bin1->sum_b - inst.b_ptr[idx1] + inst.b_ptr[idx2];
                double new_a2 = bin2->sum_a - inst.a_ptr[idx2] + inst.a_ptr[idx1];
                double new_b2 = bin2->sum_b - inst.b_ptr[idx2] + inst.b_ptr[idx1];
                double new_weight1 = new_a1 + inst.rho * sqrt(new_b1);
                double new_weight2 = new_a2 + inst.rho * sqrt(new_b2);
                if (new_weight1 <= inst.capacity && new_weight2 <= inst.capacity) {
                    double old_weight = bin1->tolWeight + bin2->tolWeight;
                    double new_weight = new_weight1 + new_weight2;
                    if (new_weight < old_weight - 1e-8) { // 改善了
                        // 执行交换
                        bin1->itemSet[idx1] = 0;
                        bin1->itemSet[idx2] = 1;
                        bin2->itemSet[idx2] = 0;
                        bin2->itemSet[idx1] = 1;
                        bin1->sum_a = new_a1;
                        bin1->sum_b = new_b1;
                        bin1->tolWeight = new_weight1;
                        bin2->sum_a = new_a2;
                        bin2->sum_b = new_b2;
                        bin2->tolWeight = new_weight2;
                        improved = true;
                        items2.erase(items2.begin() + cnt);
                        break;
                    }
                }
                ++cnt;
            }
        }
    }

    //Minimum Bin Removal
    // Step 1: find the bin with the minimum number of items
    int minIdx = -1, minItemCount = 1e10;
    for (int i = 0; i < sol.size(); ++i) {
        int count = std::accumulate(sol[i]->itemSet.begin(), sol[i]->itemSet.end(), 0);
        if (count > 0 && count < minItemCount) {
            minItemCount = count;
            minIdx = i;
        }
    }
    Bin* minBin = sol[minIdx];
    vector<int> itemsToMove;
    for (int i = 0; i < inst.n_items; ++i) {
        if (minBin->itemSet[i])
            itemsToMove.push_back(i);
    }
    // Step 2: try to insert the items into other bins
    bool removeFlag = true;
    for (int itemIdx : itemsToMove) {
        bool inserted = false;
        for (int j = 0; j < sol.size(); ++j) {
            if (j == minIdx) continue;
            Bin* target = sol[j];
            double new_a = target->sum_a + inst.a_ptr[itemIdx];
            double new_b = target->sum_b + inst.b_ptr[itemIdx];
            double new_weight = new_a + inst.rho * sqrt(new_b);
            if (new_weight <= inst.capacity + EX && target->itemSet[itemIdx] == 0) {
                // insert item
                target->itemSet[itemIdx] = 1;
                target->sum_a = new_a;
                target->sum_b = new_b;
                target->tolWeight = new_weight;
                //update information in minBin
                minBin->itemSet[itemIdx] = 0;
                minBin->sum_a -= inst.a_ptr[itemIdx];
                minBin->sum_b -= inst.b_ptr[itemIdx];
                minBin->tolWeight = minBin->sum_a + inst.rho * sqrt(minBin->sum_b);
                inserted = true;
                break;
            }
        }
        // stop when inserting failed
        if (!inserted) {
            removeFlag = false;
            break;
        }
    }
    // Step 3: if all the items are reinserted, then delete minBin
    if (removeFlag) {
        delete minBin;
        sol.erase(sol.begin() + minIdx);
    }

    ////test
    //vector<int> allItem(inst.n_items, 1);
    //for (int i = 0; i < sol.size(); ++i) {
    //    TestBinCapacity(sol[i]->itemSet, inst);
    //    for (int j = 0; j < inst.n_items; ++j) {
    //        if(sol[i]->itemSet[j] > 0)
    //            --allItem[j];
    //        if (allItem[j] < 0)
    //            cerr << "Item visit error!" << endl;
    //    }
    //}
    //int sumItem = accumulate(allItem.begin(), allItem.end(), 0);
    //if (sumItem > 0) {
    //    cerr << "Some items missed!" << endl;
    //}

    // remove the bin without items
    sol.erase(std::remove_if(sol.begin(), sol.end(), [](Bin* b) {
        bool isEmpty = std::accumulate(b->itemSet.begin(), b->itemSet.end(), 0) == 0;
        if (isEmpty) {
            delete b;  // free space
        }
        return isEmpty;
        }), sol.end());

    ////test solution correctness
    //int allItems = 0;
    //for (auto& e : sol) {
    //    allItems += std::accumulate(e->itemSet.begin(), e->itemSet.end(), 0);
    //    
    //}
    //if (allItems != inst.n_items) {
    //    cerr << "Heuristic solution wrong!" << endl;
    //}
}

//copy solutions
void CopySolutions(vector<Bin*>& preSol, vector<Bin*>& oldSol) {
    preSol.resize(oldSol.size());
    for (int i = 0; i < oldSol.size(); ++i)
        preSol[i] = new Bin(oldSol[i]);
}

//delete solution
void DelSolution(vector<Bin*>& preSol) {
    for (auto& e : preSol) {
        delete e;
        e = nullptr;
    }
}

//Initialize the population
void shuffle_instance(Instance* T) {
    if (!T || T->n_items == 0) return;

    size_t n = T->n_items;
    for (size_t i = n - 1; i > 0; --i) {
        size_t j = rand() % (i + 1); // generate a random number

        // swap a_ptr[item] nd a_ptr[j]
        double tmp = T->a_ptr[i];
        T->a_ptr[i] = T->a_ptr[j];
        T->a_ptr[j] = tmp;

        // swap b_ptr[item] and b_ptr[j]
        tmp = T->b_ptr[i];
        T->b_ptr[i] = T->b_ptr[j];
        T->b_ptr[j] = tmp;

        // swap p_ptr[item] and p_ptr[j]
        tmp = T->p_ptr[i];
        T->p_ptr[i] = T->p_ptr[j];
        T->p_ptr[j] = tmp;

        // swap p_weight[item] and p_weight[j]
        tmp = T->p_weight[i];
        T->p_weight[i] = T->p_weight[j];
        T->p_weight[j] = tmp;
    }
}
void InitializePopulation(vector<vector<Bin*>>& population, vector<Bin*>& greedyInitSol, Instance& preInstances) {
    for (int p = 0; p < g_populationSize; ++p) {
        vector<Bin*> preSol;
        if (rand() % 10 < 5) {  //get a random init solution 50%
            vector<int> insertSeq(preInstances.n_items);
            for (int i = 0; i < preInstances.n_items; ++i)
                insertSeq[i] = i;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            shuffle(insertSeq.begin(), insertSeq.end(), std::default_random_engine(seed));

            Bin* preBin = new Bin();
            preBin->itemSet.resize(preInstances.n_items);
            for (int i = 0; i < preInstances.n_items; ++i) {
                int item = insertSeq.back();
                //judge the capacity constraints
                double tmpWeight = preBin->sum_a + preInstances.a_ptr[item] +
                    preInstances.rho * sqrt(preBin->sum_b + preInstances.b_ptr[item]);

                if (tmpWeight <= preInstances.capacity + EX) {//add item to the current bin
                    preBin->itemSet[item] = 1;
                    preBin->sum_a += preInstances.a_ptr[item];
                    preBin->sum_b += preInstances.b_ptr[item];
                    insertSeq.pop_back();
                    preBin->tolWeight = tmpWeight;
                }
                else {//open a new bin
                    preSol.emplace_back(preBin);
                    preBin = new Bin();
                    preBin->itemSet.resize(preInstances.n_items);
                    --i;
                }
            }
            if (preBin->sum_a > 0)
                preSol.emplace_back(preBin);
        }
        else  //get a greedy init solution 50%
            CopySolutions(preSol, greedyInitSol);

        population.emplace_back(preSol);
    }
}


//Crossover
void Crossover(const vector<Bin*>& p1, const vector<Bin*>& p2, Instance& inst, vector<Bin*>& childSol) {
    //get two itemSets of the p1 and p2
    vector<int> items1(inst.n_items, -1), items2(inst.n_items, -1);
    int cnt = 0;
    for (auto& e : p1) {
        for (int i = 0; i < e->itemSet.size(); ++i) {
            if (e->itemSet[i]) {
                items1[cnt] = i;
                ++cnt;
            }
        }
    }
    cnt = 0;
    for (auto& e : p2) {
        for (int i = 0; i < e->itemSet.size(); ++i) {
            if (e->itemSet[i]) {
                items2[cnt] = i;
                ++cnt;
            }
        }
    }

    //random select a range from items1
    int leftIdx = rand() % items1.size();
    int rightIdx = leftIdx + rand() % (items1.size() - leftIdx);

    //extract [leftIdx, rightIdx] from items1
    vector<int> rangeSet(items1.begin() + leftIdx, items1.begin() + rightIdx);
    //vector<int> tSet1(items1.begin()+ rightIdx+1,items1.end());
    //items1.resize(leftIdx - 1);
    //items1.insert(items1.end(), tSet1.begin(), tSet1.end());

    //extract [leftIdx, rightIdx] from items2
    unordered_set<int> refSet(rangeSet.begin(), rangeSet.end());
    for (int i = 0; i < items2.size(); ++i) {
        if (refSet.find(items2[i]) != refSet.end()) {
            items2.erase(items2.begin() + i);
            --i;
        }
    }

    //reconstruct items2
    items2.insert(items2.end(), rangeSet.begin(), rangeSet.end());

    //split items in items2 greedily to form a new solution
    Bin* preBin = new Bin();
    preBin->itemSet.resize(inst.n_items);
    for (int i = 0; i < items2.size(); ++i) {
        int preItem = items2[i];
        //judge the capacity constraints
        double tmpWeight = preBin->sum_a + inst.a_ptr[preItem] +
            inst.rho * sqrt(preBin->sum_b + inst.b_ptr[preItem]);

        if (tmpWeight <= inst.capacity + EX) {//add item to the current bin
            preBin->itemSet[preItem] = 1;
            preBin->sum_a += inst.a_ptr[preItem];
            preBin->sum_b += inst.b_ptr[preItem];
            preBin->tolWeight = tmpWeight;
        }
        else {//open a new bin
            childSol.emplace_back(preBin);
            preBin = new Bin();
            preBin->itemSet.resize(inst.n_items);
            --i;
        }
    }
    if (preBin->sum_a > 0)
        childSol.emplace_back(preBin);

}

//Mutate
struct RemovedItem {
public:
    RemovedItem(double pab, int item) :p_ab(pab), itemIdx(item) {}
    RemovedItem() = default;
    bool operator < (const RemovedItem* b) const {
        return p_ab > b->p_ab; //Descending order after sorting
    }
public:
    double p_ab;
    int itemIdx;
};
void Mutate(vector<Bin*>& sol, Instance& inst) {
    int totalItems = inst.n_items;
    vector<RemovedItem*> itemToReinsert;

    // random choose some items
    for (auto& bin : sol) {
        for (int i = 0; i < totalItems; ++i) {
            if (bin->itemSet[i] && rand() % 10 < 3) { // 30% 概率移除
                itemToReinsert.push_back(
                    new RemovedItem(inst.p_ptr[i] / (inst.a_ptr[i] + inst.rho * sqrt(inst.b_ptr[i])), i)
                );
                bin->itemSet[i] = 0;
                bin->sum_a -= inst.a_ptr[i];
                bin->sum_b -= inst.b_ptr[i];
            }
        }
    }

    // remove the bin without items
    sol.erase(std::remove_if(sol.begin(), sol.end(), [](Bin* b) {
        bool isEmpty = std::accumulate(b->itemSet.begin(), b->itemSet.end(), 0) == 0;
        if (isEmpty) {
            delete b;  // free space
        }
        return isEmpty;
        }), sol.end());

    if (rand() % 10 < 2) {//20% randomly shuffle
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        shuffle(itemToReinsert.begin(), itemToReinsert.end(), std::default_random_engine(seed));
    }
    else {//80% sorted the items according p_ab
        sort(itemToReinsert.begin(), itemToReinsert.end(), [](const RemovedItem* a, const RemovedItem* b) {
            return a->p_ab > b->p_ab;  // descending order
            });
    }

    // reinsert the removed items
    for (auto& item : itemToReinsert) {
        bool inserted = false;
        multimap<double, int> feasPosRec;       //the record of the feasible insert positions, which bin to insert
        int cnt = 0;
        //record all the feasible insert positions
        for (auto& bin : sol) {
            /*double tmpWeight = bin->sum_a + inst.a_ptr[item->itemIdx] + inst.rho * sqrt(bin->sum_b + inst.b_ptr[item->itemIdx]);
            if (tmpWeight <= inst.capacity) {
                bin->itemSet[item->itemIdx] = 1;
                bin->sum_a += inst.a_ptr[item->itemIdx];
                bin->sum_b += inst.b_ptr[item->itemIdx];
                bin->tolWeight = tmpWeight;
                inserted = true;
                break;
            }*/
            double presentWeight = bin->sum_a + inst.rho * sqrt(bin->sum_b);
            double newWeight = bin->sum_a + inst.a_ptr[item->itemIdx] + inst.rho * sqrt(bin->sum_b + inst.b_ptr[item->itemIdx]);
            if (newWeight <= inst.capacity + EX)
                feasPosRec.insert({ newWeight - presentWeight, cnt });
            ++cnt;
        }
        if (!feasPosRec.empty()) {
            int insertBinIdx = -1;
            if (rand() % 10 < 5) {
                //insert the item greedily 50%
                insertBinIdx = feasPosRec.begin()->second;
            }
            else if (rand() % 10 < 6 && feasPosRec.size() > 1) {
                //insert the item using the second good position 70%
                insertBinIdx = (++feasPosRec.begin())->second;
            }
            else {
                //insert the item randomly
                int randIdx = rand() % feasPosRec.size();
                int cnt = -1;
                auto ite = feasPosRec.begin();
                while (++cnt < randIdx)
                    ++ite;
                insertBinIdx = ite->second;
            }
            sol[insertBinIdx]->itemSet[item->itemIdx] = 1;
            sol[insertBinIdx]->sum_a += inst.a_ptr[item->itemIdx];
            sol[insertBinIdx]->sum_b += inst.b_ptr[item->itemIdx];
            sol[insertBinIdx]->tolWeight =
                sol[insertBinIdx]->sum_a + inst.a_ptr[item->itemIdx] + inst.rho * sqrt(sol[insertBinIdx]->sum_b + inst.b_ptr[item->itemIdx]);;
            inserted = true;
        }

        //No place to insert the current item, then create a new bin
        if (!inserted) {
            Bin* newBin = new Bin();
            newBin->itemSet.resize(inst.n_items);
            newBin->itemSet[item->itemIdx] = 1;
            newBin->sum_a = inst.a_ptr[item->itemIdx];
            newBin->sum_b = inst.b_ptr[item->itemIdx];
            newBin->tolWeight = newBin->sum_a + inst.rho * sqrt(newBin->sum_b);
            sol.push_back(newBin);
        }
    }
    for (auto& e : itemToReinsert)
        delete e;
}

//remove solutions from population
void RomoveSolutions(vector<vector<Bin*>>& population) {
    //remove the oldest ones
    reverse(population.begin(), population.end());
    while (population.size() > g_populationSize) {
        for (auto& e : population.back())
            delete e;
        population.pop_back();
    }
}

//use hybrid genetic search to optimize the initial solution
void HGSOptimize(vector<Bin*>& bestSol, Instance& preInstances) {
    srand(static_cast<unsigned int>(time(NULL)));  // srand seed
    int outIteration = 200;
    int inIteration = 5000;
    int outCnt = 0;

    vector<vector<Bin*>> population;                //the population set
    //Initialize the population
    InitializePopulation(population, bestSol, preInstances);

    //the main loop
    while (++outCnt <= outIteration) {
        vector<Bin*> preSol;
        //50% use bestSol to optimize
        if (rand() % 10 < 5) {
            DelSolution(preSol);
            CopySolutions(preSol, bestSol);
        }
        else {
            //randomly choose two parent solutions
            size_t pIdx1 = rand() % population.size();
            size_t pIdx2 = rand() % (population.size() - 1);
            if (pIdx2 >= pIdx1)
                ++pIdx2; // idx1 differs from idx2
            /*use the genetic algorithm */
            //Crossover
            Crossover(population[pIdx1], population[pIdx2], preInstances, preSol);
        }
        if (preSol.size() < bestSol.size()) {
            cout << "Crossover gets a new heuristic solution with " << preSol.size() << " bins" << endl;
            DelSolution(bestSol);
            CopySolutions(bestSol, preSol);
            population.emplace_back(preSol);
            for (int i = 0; i < 10; ++i) {
                vector<Bin*> tmpSol;
                CopySolutions(tmpSol, preSol);
                population.emplace_back(tmpSol);
            }
            continue;
        }

        //mutate
        Mutate(preSol, preInstances);
        if (preSol.size() < bestSol.size()) {
            ////test the capacity
            //for (auto& e : preSol) {
            //    if (e->tolWeight > g_instance.capacity) {
            //        int a = 0;
            //    }
            //}
            cout << "Mutate gets a new heuristic solution with " << preSol.size() << " bins" << endl;
            DelSolution(bestSol);
            CopySolutions(bestSol, preSol);
            population.emplace_back(preSol);
            for (int i = 0; i < 10; ++i) {
                vector<Bin*> tmpSol;
                CopySolutions(tmpSol, preSol);
                population.emplace_back(tmpSol);
            }
            continue;
        }

        int inCnt = 0;
        bool improveFlag = false;
        while (++inCnt <= inIteration) {
            /*use the local search*/
            LocalSearch(preSol, preInstances);
            if (preSol.size() < bestSol.size()) {
                ////test the capacity
                //for (auto& e : preSol) {
                //    if (e->tolWeight > g_instance.capacity) {
                //        int a = 0;
                //    }
                //}
                cout << "LS gets a new heuristic solution with " << preSol.size() << " bins" << endl;
                DelSolution(bestSol);
                CopySolutions(bestSol, preSol);
                improveFlag = true;
                population.emplace_back(preSol);
                for (int i = 0; i < 10; ++i) {
                    vector<Bin*> tmpSol;
                    CopySolutions(tmpSol, preSol);
                    population.emplace_back(tmpSol);
                }
                break;
            }
        }
        //free preSol
        if (!improveFlag)
            DelSolution(preSol);

        //remove solutions from population
        if (population.size() >= g_populationMaxSize)
            RomoveSolutions(population);
    }

    //free population
    for (auto& e : population) {
        for (auto& e1 : e)
            delete e1;
    }
}

//call the heuristic to solve the submodular bin packing problem
int HeuristicForBinpacking(vector<Bin*>& bestSol) {
    auto startTime = chrono::high_resolution_clock::now();
    Instance preInstances = { 0 };
    CopyInstance(preInstances, g_instance);
    preInstances.n_items = g_instance.n_items;
    preInstances.rho = g_instance.rho;
    preInstances.capacity = g_instance.capacity;

    // sort the items according to the p/(a+rho*sqrt(b))
    ItemIndex* indicies = nullptr;
    sort_instance_by_p_ab(preInstances, indicies);

    //get initial solution
    GenerateInitialSol(bestSol, preInstances);

    //use HGS to optimize
    HGSOptimize(bestSol, preInstances);


    //adjust the indicies according to indicies
    for (auto& e : bestSol) {
        vector<int> preItemSet(preInstances.n_items);
        for (int i = 0; i < preInstances.n_items; ++i) {
            if (e->itemSet[i] == 1)
                preItemSet[indicies[i].index] = 1;
        }
        e->itemSet = preItemSet;
    }
    free(indicies);
    free(preInstances.a_ptr);
    free(preInstances.b_ptr);
    free(preInstances.p_ptr);
    free(preInstances.p_weight);
    preInstances.a_ptr = preInstances.b_ptr = preInstances.p_ptr = preInstances.p_weight = nullptr;

    auto endTime = chrono::high_resolution_clock::now();
    g_heuristicTime += chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0;
    //std::this_thread::sleep_for(std::chrono::seconds(2)); // 睡眠 2 秒
    return bestSol.size();
}