# @ Time : 2020/12/15
# @ Author: Tu Yanlun
# @ Team mates: Xie Yuzhang, Yu Jingwei

import numpy as np
from ReadData import get_data, job_task_relationship, get_index
from ProcessData import Mapper, D_kis, C_kij, E_kij, get_M, A_j, W_kij, Precedence
from LPsolver import LPsolver, compute_max_phi, fix_x_kij, update_resource_capacity

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
    Precedence_instance = Precedence(Task_List, mapper_instance, Job_Precedence=Job_Precedence)
    precedence = Precedence_instance.get_precedence()
    W_kij_instance = W_kij(c_kij=C_kij_instance.get_c_kij(), e_kij=E_kij_instance.get_e_kij(), precedence=precedence)

    max_k, max_i, max_j = C_kij_instance.get_c_kij().shape

    fix_kij_tuple_list = []
    visited_k = set()
    init_x_kij = np.zeros((max_k, max_i, max_j))
    x_kij = None
    while (len(visited_k) < max_k-2):
        print("the {} iteration".format(len(visited_k)))
        print("------------------------------------------------------")
        W_kij_instance.compute_w_kij()
        w_kij = W_kij_instance.get_w_kij()

        x_kij = LPsolver(c_kij=C_kij_instance.get_c_kij(), e_kij=E_kij_instance.get_e_kij(), M=M,
                         a_j=A_j_instance.get_a_j(),
                         job_task_mapping=mapper_instance.get_job_task_idx_mapping(),
                         fix_kij_tuple_list=fix_kij_tuple_list,
                         init_x_kij=init_x_kij,
                         w_kij=w_kij,
                         precedence=precedence,
                         b_sj=C_kij_instance.get_b_sj())

        res_k, res_i, res_j = compute_max_phi(x_kij=x_kij, c_kij=C_kij_instance.get_c_kij(),
                                              e_kij=E_kij_instance.get_e_kij(), visited_k=visited_k)

        visited_k.add(res_k)

        # update w_kij
        W_kij_instance.update_x_kij(x_kij)

        fix_x_kij(fix_kij_tuple_list=fix_kij_tuple_list, res_k=res_k, res_i=res_i, max_j=max_j, init_x_kij=init_x_kij,
                  x_kij=x_kij)

        update_resource_capacity(res_j=res_j, A_j=A_j_instance)

    final_x_kij = np.zeros((max_k, max_i, max_j))
    for k in range(max_k):
        for i in range(max_i):
            for j in range(max_j):
                if (k, i, j) in fix_kij_tuple_list:
                    final_x_kij[k][i][j] = init_x_kij[k][i][j]
                else:
                    final_x_kij[k][i][j] = x_kij[k][i][j]
    np.save("final_x.npy", final_x_kij)
    print("done.....................................................")
