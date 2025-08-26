#include "compactModel.h"
#include "labelSetting.h"

using namespace std;

// call the gurobi solver to solve the compact SMBP model
int SolveCompactSMBPModel(Args& args, double binLB, int binUB) {
    //read instance
    Instance instance = { 0 };
    if (instance_parse(&instance, args.input_file)) {
        cerr << "Error reading instance file" << endl;
        return -1;
    }
    if (instance.rho < 1 && !args.direct_rho)
        instance.rho = sqrt(instance.rho) / sqrt(1 - instance.rho);
    ofstream outPut(args.output_file, ios::app);

    //feasible flag, 0 indicates infeasible, 1 indicates integer optimal solution, 2 indicates feasible solution but not optimal
    int feasFlag = true;
    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        GRBModel model = GRBModel(env);
        //variable v specify if the item i is inserted into the knapsack j
        GRBVar** v = new GRBVar * [instance.n_items];
        for (int i = 0; i < instance.n_items; ++i)
            v[i] = model.addVars(binUB, GRB_BINARY);
        //variable y specify if the bin j is used
        GRBVar* y = model.addVars(binUB, GRB_BINARY);
        //variable t is the auxiliary variable for each SOCP constraints
        GRBVar* t = model.addVars(binUB, GRB_CONTINUOUS);


        //set objective function
        GRBLinExpr obj = 0;
        for (int i = 0; i < binUB; ++i)
            obj += y[i];
        model.setObjective(obj, GRB_MINIMIZE);

        //set SOCP constraints
        for (int j = 0; j < binUB; ++j) {//for each bin
            GRBLinExpr lhs = 0;
            for (int i = 0; i < instance.n_items; ++i)
                lhs += instance.a_ptr[i] * v[i][j]; //a_i * v_ij, sum a part
            GRBLinExpr rhs = 0;
            for (int i = 0; i < instance.n_items; ++i)
                rhs += instance.b_ptr[i] * v[i][j]; //b_j * v_ij, sum b part
            model.addQConstr(t[j] * t[j] == rhs, "t_squared_constraint");

            model.addConstr(lhs + instance.rho * t[j] <= instance.capacity * y[j], "SOCP_constraint");
        }

        //each item can only be inserted into one bin
        for (int i = 0; i < instance.n_items; ++i) {
            GRBLinExpr itemSum = 0;
            for (int j = 0; j < binUB; ++j)
                itemSum += v[i][j];
            model.addConstr(itemSum == 1, "item_constraint_" + to_string(i)); //each item must be in exactly one bin
        }

        //test, set the lower and upper bound for the bin number
        model.addConstr(obj >= binLB);
        model.addConstr(obj <= binUB);

        // set the model precision parameters
        model.set(GRB_DoubleParam_MIPGap, 1e-6);           
        model.set(GRB_DoubleParam_FeasibilityTol, EX);     
        model.set(GRB_DoubleParam_OptimalityTol, EX);      
        model.set(GRB_DoubleParam_IntFeasTol, EX);         
        model.set(GRB_IntParam_Presolve, 2);               
        model.set(GRB_IntParam_Method, 2);                 

        model.set(GRB_DoubleParam_TimeLimit, 60);       //set maximum testing time is 60s

        // solve the model
        model.optimize();

        cout << "***************The results obtained by Gurobi**************" << endl;
        outPut << "\n\n***************The results obtained by Gurobi**************" << endl;
        // output the results
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            feasFlag = 1;
            g_binpackingUB = model.get(GRB_DoubleAttr_ObjVal);
        }
        else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {//feasible but not optimal
            cout << "Reach time limitation: " << endl;
            outPut << "Reach time limitation: " << endl;
            feasFlag = 2;
            if (model.get(GRB_DoubleAttr_ObjVal) < g_binpackingUB) {
                g_binpackingUB = model.get(GRB_DoubleAttr_ObjVal);//update upper bound
            }
        }
        else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
            cout << "Model infeasible!" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            outPut << "Model infeasible!" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            feasFlag = 0;
            return feasFlag;
        }
        cout << "The computational time is: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        cout << "The number of bins used is (upper bound): " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        cout << "Best Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
        cout << "MIP Gap: " << model.get(GRB_DoubleAttr_MIPGap) << endl;
        outPut << "The computational time is: " << model.get(GRB_DoubleAttr_Runtime) << endl;
        outPut << "The number of bins used is (upper bound): " << model.get(GRB_DoubleAttr_ObjVal) << endl;
        outPut << "Best Bound: " << model.get(GRB_DoubleAttr_ObjBound) << endl;
        outPut << "MIP Gap: " << model.get(GRB_DoubleAttr_MIPGap) << endl;
        //output the capacity of each bin
        for (int j = 0; j < binUB; ++j) {
            if (y[j].get(GRB_DoubleAttr_X) >= 0.99) {
                double lhs = 0;
                for (int i = 0; i < instance.n_items; ++i)
                    lhs += instance.a_ptr[i] * v[i][j].get(GRB_DoubleAttr_X); //a_i * v_ij, sum a part
                double rhs = 0;
                for (int i = 0; i < instance.n_items; ++i)
                    rhs += instance.b_ptr[i] * v[i][j].get(GRB_DoubleAttr_X); //b_j * v_ij, sum b part
                cout << "The bin use capacity " << lhs + instance.rho * sqrt(rhs) << endl;
                outPut << "The bin use capacity " << lhs + instance.rho * sqrt(rhs) << endl;
            }
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi error code: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        outPut << "Gurobi error code: " << e.getErrorCode() << std::endl;
        outPut << e.getMessage() << std::endl;
        feasFlag = 2;
        return feasFlag;
    }
    catch (exception& e) {
        std::cout << "Other errors:" << e.what() << std::endl;
        outPut << "Other errors:" << e.what() << std::endl;
        feasFlag = 2;
        return feasFlag;
    }

    outPut.close();
    //free space
    instance_free(&instance);
    return feasFlag;
}

//call the gurobi solver to solve the compact SOCP model of submodular knapsack problem
void SolveCompactKnapsackModel(Args& args) {
    //read instance
    Instance instance = { 0 };
    if (instance_parse(&instance, args.input_file)) {
        cerr << "Error reading instance file" << endl;
        return;
    }
    if (instance.rho < 1 && !args.direct_rho)
        instance.rho = sqrt(instance.rho) / sqrt(1 - instance.rho);
    ofstream outPut(args.output_file, ios::app);

    try {
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "./gurobi.log");
        env.start();
        GRBModel model = GRBModel(env);

        //variable specify if the item i is selected
        GRBVar* x = model.addVars(instance.n_items, GRB_BINARY);
        GRBVar t = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        //set objective function
        GRBLinExpr obj = 0;
        for (int i = 0; i < instance.n_items; ++i)
            obj += instance.p_ptr[i] * x[i];
        model.setObjective(obj, GRB_MAXIMIZE);  

        //set SOCP constraints
        GRBLinExpr lhs = 0;
        for (int j = 0; j < instance.n_items; ++j)
            lhs += instance.a_ptr[j] * x[j];
        model.addConstr(lhs + instance.rho * t <= instance.capacity, "SOCP_constraint");
        GRBLinExpr rhs = 0;
        for (int j = 0; j < instance.n_items; ++j)
            rhs += instance.b_ptr[j] * x[j];
        model.addQConstr(t * t == rhs, "t_squared_constraint");

        // set precision
        model.set(GRB_DoubleParam_MIPGap, 1e-6);           
        model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);   
        model.set(GRB_DoubleParam_OptimalityTol, 1e-9);    
        model.set(GRB_DoubleParam_IntFeasTol, 1e-9);       
        model.set(GRB_IntParam_Presolve, 2);               
        model.set(GRB_IntParam_Method, 2);                 

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);        //set maximum solution time

        // solve the model
        model.optimize();

        cout << "***************The results obtained by Gurobi**************" << endl;
        outPut << "\n\n***************The results obtained by Gurobi**************" << endl;
        // output the results
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            cout << "The optimal objective value is :" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            cout << "The solution time is :" << model.get(GRB_DoubleAttr_Runtime) << endl;
            cout << "The chosen items are as follows: " << endl;
            outPut << "The optimal objective value is :" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            outPut << "The solution time is :" << model.get(GRB_DoubleAttr_Runtime) << endl;
            outPut << "The chosen items are as follows: " << endl;

            int itemNum = 0;
            for (int i = 0; i < instance.n_items; i++) {
                if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
                    ++itemNum;
                    cout << i << ",";
                    outPut << i << ",";
                }
            }
            cout << endl;
            cout << "The total number of chosen items is: " << itemNum << endl;
            outPut << endl;
            outPut << "The total number of chosen items is: " << itemNum << endl;
            std::cout << "t = " << t.get(GRB_DoubleAttr_X) << std::endl;
        }
        else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            cout << "Time limit reached. Returning the best solution found so far." << endl;
            outPut << "Time limit reached. Returning the best solution found so far." << endl;

            // Retrieve the best solution found up to the timeout
            double bestObjVal = model.get(GRB_DoubleAttr_ObjVal);
            cout << "The best objective value found is: " << bestObjVal << endl;
            outPut << "The best objective value found is: " << bestObjVal << endl;

            int itemNum = 0;
            for (int i = 0; i < instance.n_items; i++) {
                if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
                    ++itemNum;
                    cout << i << ",";
                    outPut << i << ",";
                }
            }
            cout << endl;
            cout << "The total number of chosen items is: " << itemNum << endl;
            outPut << endl;
            outPut << "The total number of chosen items is: " << itemNum << endl;
            std::cout << "t = " << t.get(GRB_DoubleAttr_X) << std::endl;
        }
        else {
            std::cout << "No optimal solution found" << std::endl;
            outPut << "No optimal solution found" << std::endl;
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi error code: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        outPut << "Gurobi error code: " << e.getErrorCode() << std::endl;
        outPut << e.getMessage() << std::endl;
        return;
    }
    catch (exception& e) {
        std::cout << "Other errors:" << e.what() << std::endl;
        outPut << "Other errors:" << e.what() << std::endl;
        return;
    }

    outPut.close();
    //free space
    instance_free(&instance);
}


//call the gurobi solver to solve the pricing problem
int SolveCompactPricingModel(
    vector<Bin*>& newCols,
    ofstream& outPut,
    vector<Stage>& preStages,
    double& preTime,
    double& binNumDual,
    vector<double>& preciseDuals
) {
    auto startTime = chrono::high_resolution_clock::now();
    //sort the instance according to p/(a+sqrt(b))
    ItemIndex* indicesRec = nullptr;
    sort_instance_by_p_ab(g_instance, indicesRec);

    /*adjust instances according to the preStages*/
    unordered_map<int, ItemIndex> removedIndicesRec;
    AdjustInstancesPerStages(preStages, indicesRec, removedIndicesRec);

    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        env.set(GRB_IntParam_OutputFlag, 0); // trun off the output
        GRBModel model = GRBModel(env);

        //variable specify if the item i is selected
        GRBVar* x = model.addVars(g_instance.n_items, GRB_BINARY);
        GRBVar t = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        //set objective function
        GRBLinExpr obj = 0;
        for (int i = 0; i < g_instance.n_items; ++i)
            obj += g_instance.p_ptr[i] * x[i];
        model.setObjective(obj, GRB_MAXIMIZE);  

        //set SOCP constraints
        GRBLinExpr lhs = 0;
        for (int j = 0; j < g_instance.n_items; ++j)
            lhs += g_instance.a_ptr[j] * x[j];
        model.addConstr(lhs + g_instance.rho * t <= g_instance.capacity, "SOCP_constraint");
        GRBLinExpr rhs = 0;
        for (int j = 0; j < g_instance.n_items; ++j)
            rhs += g_instance.b_ptr[j] * x[j];
        model.addQConstr(t * t == rhs, "t_squared_constraint");

        // set precision
        model.set(GRB_DoubleParam_MIPGap, 1e-6);         
        model.set(GRB_DoubleParam_FeasibilityTol, EX);   
        model.set(GRB_DoubleParam_OptimalityTol, EX);    
        model.set(GRB_DoubleParam_IntFeasTol, EX);       
        model.set(GRB_IntParam_Presolve, 2);             
        model.set(GRB_IntParam_Method, 2);               

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);         

        // solve the model
        model.optimize();

        // obtain the best results
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            double objVal = model.get(GRB_DoubleAttr_ObjVal);
            int solCount = model.get(GRB_IntAttr_SolCount);
            cout << "The optimal objective value is :" << objVal << endl;
            cout << "Total number of solution found is :" << solCount << endl;

            for (int i = 0; i < solCount; ++i) {
                model.set(GRB_IntParam_SolutionNumber, i);  // obtain the present solution
                objVal = model.get(GRB_DoubleAttr_ObjVal);
                if (1.0 - (objVal + binNumDual) < -1e-6) {
                    newCols.push_back(new Bin());
                    newCols.back()->itemSet.resize(g_tolItemNum);
                    newCols.back()->tolProfit = objVal;
                    for (int i = 0; i < g_instance.n_items; ++i) {
                        if (x[i].get(GRB_DoubleAttr_X) >= 0.99)
                            newCols.back()->itemSet[indicesRec[i].index] = 1;
                    }
                }
            }

        }
        else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            cout << "Time limit reached. Returning the best solution found so far." << endl;
            outPut << "Time limit reached. Returning the best solution found so far." << endl;

            // Retrieve the best solution found up to the timeout
            double bestObjVal = model.get(GRB_DoubleAttr_ObjVal);
            cout << "The best objective value found is: " << bestObjVal << endl;
            outPut << "The best objective value found is: " << bestObjVal << endl;

            int itemNum = 0;
            for (int i = 0; i < g_instance.n_items; i++) {
                if (x[i].get(GRB_DoubleAttr_X) >= 0.9) {
                    ++itemNum;
                    cout << i << ",";
                    outPut << i << ",";
                }
            }
            cout << endl;
            cout << "The total number of chosen items is: " << itemNum << endl;
            outPut << endl;
            outPut << "The total number of chosen items is: " << itemNum << endl;
            std::cout << "t = " << t.get(GRB_DoubleAttr_X) << std::endl;
        }
        else {
            std::cout << "No optimal solution found" << std::endl;
            outPut << "No optimal solution found" << std::endl;
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi error code: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        outPut << "Gurobi error code: " << e.getErrorCode() << std::endl;
        outPut << e.getMessage() << std::endl;
        return 0;
    }
    catch (exception& e) {
        std::cout << "Other errors:" << e.what() << std::endl;
        outPut << "Other errors:" << e.what() << std::endl;
        return 0;
    }

    outPut.close();
    //free space
    instance_free(&g_instance);
    g_instance.a_ptr = g_instance.b_ptr = g_instance.p_ptr = g_instance.p_weight = nullptr;
    return newCols.size();
}