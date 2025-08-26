#include "compactModel.h"

using namespace std;

//call the gurobi solver to solve the linear relaxation of the knapsack problem
void SolveLinearKnapsack(Instance& instance, int itemIdx, double& preWeight, double& dualProfit) {
    try {
        GRBEnv env = GRBEnv(true);
        env.start();
        env.set(GRB_IntParam_OutputFlag, 0);  // 禁用所有输出
        GRBModel model = GRBModel(env);

        //variable specify if the item i is selected
        GRBVar* x = model.addVars(itemIdx, GRB_CONTINUOUS);
        GRBVar t = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);

        //set objective function
        GRBLinExpr obj = 0;
        for (int i = 0; i < itemIdx; ++i)
            obj += instance.p_ptr[i] * x[i];
        model.setObjective(obj, GRB_MAXIMIZE);  // 设置目标函数为最大化

        //set SOCP constraints
        GRBLinExpr lhs = 0;
        for (int j = 0; j < itemIdx; ++j)
            lhs += instance.a_ptr[j] * x[j];
        model.addConstr(lhs <= preWeight, "capacity_constraint");

        model.set(GRB_DoubleParam_TimeLimit, 7200.0);         // 设置最大求解时间为60秒

        // solve the model
        model.optimize();

        // output the results
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            //cout << "The optimal objective value is :" << model.get(GRB_DoubleAttr_ObjVal) << endl;
            //cout << "The solution time is :" << model.get(GRB_DoubleAttr_Runtime) << endl;

            dualProfit = model.get(GRB_DoubleAttr_ObjVal);
        }
        else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            cout << "Time limit reached for linear relaxation of the knapsack." << endl;

            // Retrieve the best solution found up to the timeout
            dualProfit = model.get(GRB_DoubleAttr_ObjVal);
        }
        else {
            std::cout << "Linear relaxation of the knapsack: no optimal solution found" << std::endl;
        }
    }
    catch (GRBException& e) {
        std::cout << "Gurobi error code: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        return;
    }
    catch (exception& e) {
        std::cout << "Other errors:" << e.what() << std::endl;
        return;
    }
    //cout << "The best objective value found is: " << dualProfit << endl;
}