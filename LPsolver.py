import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB
from ReadData import get_data, job_task_relationship, get_index
import math
from ProcessData import Mapper, D_kis, C_kij, E_kij, get_M, A_j


def LPsolver(c_kij, e_kij, M, a_j, job_task_mapping):
    model = gp.Model("LPsolver")

    max_k, max_i, max_j = c_kij.shape
    kij_tuple_list = []
    for k in range(max_k):
        for i in range(max_i):
            for j in range(max_j):
                kij_tuple_list.append((k, i, j))
    ki_tuple_list = []
    for k in range(max_k):
        for i in range(max_i):
            ki_tuple_list.append((k, i))
    # create a 3-D array of binary variables
    # lambda [k][i][j]
    la_0 = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="la_0")
    la_1 = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="la_1")
    x = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="x")
    model.update()

    # print(gp.quicksum(la_0[k, i, j] + math.pow(M, c_kij[k][i][j] + e_kij[k][i][j]) * la_1[k, i, j] for k, i, j in kij_tuple_list))
    # 需要优化的目标函数
    model.setObjective(gp.quicksum(
        la_0[k, i, j] + math.pow(M, c_kij[k][i][j] + e_kij[k][i][j]) * la_1[k, i, j] for k, i, j in kij_tuple_list),
        GRB.MINIMIZE)

    # constraints of paper
    # constraint 11-1
    model.addConstrs(x[k, i, j] == la_1[k, i, j] for k, i, j in kij_tuple_list)
    # constraint 11-2
    model.addConstrs(la_0[k, i, j] + la_1[k, i, j] == 1 for k, i, j in kij_tuple_list)
    # constraint 11-3 already satisfied
    # constraint 5
    for j in range(max_j):
        model.addConstr(gp.quicksum(x[k, i, j] for k, i in ki_tuple_list) <= a_j[j])
    # constraint 6
    for k in range(max_k):
        for i in range(max_i):
            if i in job_task_mapping[k]:
                model.addConstr(gp.quicksum(x[k, i, j] for j in range(max_j)) == 1)

    # additional constraint due to my design of data structure
    for k in range(max_k):
        for i in range(max_i):
            for j in range(max_j):
                if i not in job_task_mapping[k]:
                    model.addConstr(x[k, i, j] == 0)

    model.optimize()

    result_x = np.zeros((max_k,max_i,max_j))
    for k,i,j in kij_tuple_list:
        result_x[k][i][j] = x[k,i,j].x
    return result_x
    # following two sentences are used to debug model infeasible
    # model.computeIIS()
    # model.write("model1.ilp")

def compute_fi():
    pass


if __name__ == "__main__":
    Inter_Link, Job_Data, Execution_Time, Job_Precedence, Slot_Number, Data_Partition = get_data()
    Job_List, Task_List, Job_Task_List, Job_Task_Dict = job_task_relationship()
    get_index(Job_Data, Job_Task_List)

    mapper_instance = Mapper(Job_List, Task_List, Slot_Number, Data_Partition, Job_Task_Dict)
    D_kis_instance = D_kis(mapper_instance, Job_Data)
    C_kij_instance = C_kij(D_kis_instance, mapper_instance, Inter_Link)
    C_kij_instance.compute_b_sj()
    C_kij_instance.compute_c_kij()
    E_kij_instance = E_kij(Execution_Time, C_kij_instance.get_c_kij(), Job_Task_List, mapper_instance)
    M = get_M(mapper_instance, Job_List, Job_Task_Dict)
    A_j_instance = A_j(Slot_Number)
    LPsolver(c_kij=C_kij_instance.get_c_kij(), e_kij=E_kij_instance.get_e_kij(), M=M, a_j=A_j_instance.get_a_j(),
             job_task_mapping=mapper_instance.get_job_task_idx_mapping())

    print("pause")
