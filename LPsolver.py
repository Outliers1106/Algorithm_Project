# @ Time : 2020/12/15
# @ Author: Tu Yanlun
# @ Team mates: Xie Yuzhang, Yu Jingwei

import numpy as np
import gurobipy as gp
from gurobipy import GRB
import math



def LPsolver(c_kij, e_kij, M, a_j, job_task_mapping, fix_kij_tuple_list, init_x_kij):
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
    # lambda [k][i][j] x 2
    la_0 = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="la_0")
    la_1 = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="la_1")
    x = model.addVars(max_k, max_i, max_j, vtype=GRB.BINARY, name="x")
    model.update()

    # print(gp.quicksum(la_0[k, i, j] + math.pow(M, c_kij[k][i][j] + e_kij[k][i][j]) * la_1[k, i, j] for k, i, j in kij_tuple_list))
    # 需要优化的目标函数
    model.setObjective(gp.quicksum(
        la_0[k, i, j] + math.pow(M, c_kij[k][i][j] + e_kij[k][i][j]) * la_1[k, i, j] for k, i, j in kij_tuple_list if
        (k, i, j) not in fix_kij_tuple_list),
        GRB.MINIMIZE)

    # constraints of paper
    # constraint 11-1
    model.addConstrs(x[k, i, j] == la_1[k, i, j] for k, i, j in kij_tuple_list if (k, i, j) not in fix_kij_tuple_list)
    # constraint 11-2
    model.addConstrs(
        la_0[k, i, j] + la_1[k, i, j] == 1 for k, i, j in kij_tuple_list if (k, i, j) not in fix_kij_tuple_list)
    # constraint 11-3 already satisfied
    # constraint 5
    for j in range(max_j):
        model.addConstr(
            gp.quicksum(x[k, i, j] for k, i in ki_tuple_list if (k, i, j) not in fix_kij_tuple_list) <= a_j[j])
    # constraint 6
    for k in range(max_k):
        for i in range(max_i):
            if i in job_task_mapping[k]:
                model.addConstr(
                    gp.quicksum(x[k, i, j] for j in range(max_j) ) == 1)

    # additional constraint due to my design of data structure
    for k in range(max_k):
        for i in range(max_i):
            for j in range(max_j):
                if i not in job_task_mapping[k] and (k, i, j) not in fix_kij_tuple_list:
                    model.addConstr(x[k, i, j] == 0)

    # set fixed x kij (optimization function does not include them)
    for k, i, j in fix_kij_tuple_list:
        model.addConstr(x[k, i, j] == init_x_kij[k][i][j])

    model.optimize()

    result_x = np.zeros((max_k, max_i, max_j))

    for k, i, j in kij_tuple_list:
        result_x[k][i][j] = x[k, i, j].x
    return result_x
    # following two sentences are used to debug model infeasible
    # model.computeIIS()
    # model.write("model1.ilp")


# previously computed job k on phi is not computed anymore
def compute_max_phi(x_kij, c_kij, e_kij, visited_k):
    max_k, max_i, max_j = np.shape(x_kij)
    kij_tuple_list = []
    for k in range(max_k):
        for i in range(max_i):
            for j in range(max_j):
                kij_tuple_list.append((k, i, j))
    max_phi = -1
    res_k, res_i, res_j = -1, -1, -1
    for k, i, j in kij_tuple_list:
        # a trick : do not compute phi on previous computed job k, then we can get completion time by another job
        if k in visited_k:
            continue
        cur_phi = x_kij[k][i][j] * (c_kij[k][i][j] + e_kij[k][i][j])
        if cur_phi > max_phi:
            max_phi = cur_phi
            res_k, res_i, res_j = k, i, j

    return res_k, res_i, res_j


def fix_x_kij(fix_kij_tuple_list, res_k, res_i, max_j, x_kij, init_x_kij):
    for j in range(max_j):
        fix_kij_tuple_list.append((res_k, res_i, j))
        init_x_kij[res_k][res_i][j] = x_kij[res_k][res_i][j]


def update_resource_capacity(res_j, A_j):
    A_j.decrease_a_j(res_j)



# if __name__ == "__main__":
#     Inter_Link, Job_Data, Execution_Time, Job_Precedence, Slot_Number, Data_Partition = get_data()
#     Job_List, Task_List, Job_Task_List, Job_Task_Dict = job_task_relationship()
#     get_index(Job_Data, Job_Task_List)
#
#     mapper_instance = Mapper(Job_List, Task_List, Slot_Number, Data_Partition, Job_Task_Dict)
#     D_kis_instance = D_kis(mapper_instance, Job_Data)
#     C_kij_instance = C_kij(D_kis_instance, mapper_instance, Inter_Link)
#     C_kij_instance.compute_b_sj()
#     C_kij_instance.compute_c_kij()
#     E_kij_instance = E_kij(Execution_Time, C_kij_instance.get_c_kij(), Job_Task_List, mapper_instance)
#     M = get_M(mapper_instance, Job_List, Job_Task_Dict)
#     A_j_instance = A_j(Slot_Number)
#
#     max_k, max_i, max_j = C_kij_instance.get_c_kij().shape
#
#     fix_kij_tuple_list = []
#     visited_k = set()
#     init_x_kij = np.zeros((max_k, max_i, max_j))
#     x_kij = None
#     while (len(visited_k) < max_k):
#         print("the {} iteration".format(len(visited_k)))
#         print("------------------------------------------------------")
#         x_kij = LPsolver(c_kij=C_kij_instance.get_c_kij(), e_kij=E_kij_instance.get_e_kij(), M=M,
#                          a_j=A_j_instance.get_a_j(),
#                          job_task_mapping=mapper_instance.get_job_task_idx_mapping(),
#                          fix_kij_tuple_list=fix_kij_tuple_list,
#                          init_x_kij=init_x_kij)
#
#         res_k, res_i, res_j = compute_max_phi(x_kij=x_kij, c_kij=C_kij_instance.get_c_kij(),
#                                           e_kij=E_kij_instance.get_e_kij(), visited_k=visited_k)
#
#         visited_k.add(res_k)
#
#         fix_x_kij(fix_kij_tuple_list=fix_kij_tuple_list, res_k=res_k, res_i=res_i, max_j=max_j, init_x_kij=init_x_kij,
#                   x_kij=x_kij)
#
#         update_resource_capacity(res_j=res_j, A_j=A_j_instance)
#
#     final_x_kij = np.zeros((max_k, max_i, max_j))
#     for k in range(max_k):
#         for i in range(max_i):
#             for j in range(max_j):
#                 if (k,i,j) in fix_kij_tuple_list:
#                     final_x_kij[k][i][j] = init_x_kij[k][i][j]
#                 else:
#                     final_x_kij[k][i][j] = x_kij[k][i][j]
#     np.save("final_x.npy",final_x_kij)
#     print("done.....................................................")
