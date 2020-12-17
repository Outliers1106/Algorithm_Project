# @ Time : 2020/12/15
# @ Author: Tu Yanlun
# @ Team mates: Xie Yuzhang, Yu Jingwei

import numpy as np
from utils import TwoWayDict
import pandas as pd


# from ReadData import get_data, job_task_relationship, get_index
# import math


class Mapper:
    def __init__(self, Job_List, Task_List, Slot_Number, Data_Partition, Job_Task_Dict):
        self.job = TwoWayDict()
        self.task = TwoWayDict()
        self.datacenter = TwoWayDict()
        self.slots = dict()  # 每个datacenter有几个slot
        self.data_partition2datacenter = dict()
        # self.datacenter2data_partition = dict()
        for i in range(len(Job_List)):
            self.job[Job_List[i]] = i
        for i in range(len(Task_List)):
            self.task[Task_List[i]] = i
        for i in range(len(Slot_Number['DC'])):
            self.datacenter[Slot_Number['DC'][i]] = i
        for i in range(len(Slot_Number['Num of Slots'])):
            self.slots[i] = Slot_Number['Num of Slots'][i]  # 例，0对应2
            self.slots[self.datacenter[i]] = Slot_Number['Num of Slots'][i]  # 例，DC1对应2
        for i in range(len(Data_Partition['Data Partition'])):
            data_partition_name = Data_Partition['Data Partition'][i]
            datacenter_name = Data_Partition['Location'][i]
            self.data_partition2datacenter[data_partition_name] = datacenter_name
            # self.datacenter2data_partition[datacenter_name] = data_partition_name

        self.job_list = Job_List
        self.task_list = Task_List
        self.datacenter_list = list()
        for i in range(len(Slot_Number['DC'])):
            self.datacenter_list.append(Slot_Number['DC'][i])
        self.data_partition_list = list()
        for i in range(len(Data_Partition['Data Partition'])):
            self.data_partition_list.append(Data_Partition['Data Partition'][i])

        self.job_task_idx_mapping = {}
        for i in range(len(Job_List)):
            job_name = Job_List[i]
            job_id = self.get_job(job_name)
            self.job_task_idx_mapping[job_id] = list()
            for task_name in Job_Task_Dict[job_name]:
                task_id = self.get_task(task_name)
                self.job_task_idx_mapping[job_id].append(task_id)

    def get_job(self, key):
        return self.job[key]

    def get_task(self, key):
        return self.task[key]

    def get_datacenter(self, key):
        return self.datacenter[key]

    def get_slots(self, key):
        return self.slots[key]

    def get_data_partition2data_center(self, key):
        return self.data_partition2datacenter[key]

    def get_job_list(self):
        return self.job_list

    def get_task_list(self):
        return self.task_list

    def get_datacenter_list(self):
        return self.datacenter_list

    def get_data_partition_list(self):
        return self.data_partition_list

    def get_job_task_idx_mapping(self):
        # a dict
        return self.job_task_idx_mapping


class D_kis:  # checked: OK
    def __init__(self, mapper, Job_Data):
        self.mapper = mapper
        job_list = mapper.get_job_list()
        task_list = mapper.get_task_list()
        datacenter_list = mapper.get_datacenter_list()
        data_partition_list = mapper.get_data_partition_list()
        max_k = len(job_list)
        max_i = len(task_list)
        # max_s = len(data_partition_list)
        max_s = len(datacenter_list)

        max_datapartition = len(data_partition_list)

        # assign task_id to d[k][i][s]
        self.d_kis = np.zeros((max_k, max_i, max_s))

        for ss in range(max_datapartition):
            data_partition_name = data_partition_list[ss]
            datacenter_name = mapper.get_data_partition2data_center(data_partition_name)
            s = mapper.get_datacenter(datacenter_name)
            for i in range(max_i):
                job_name = Job_Data.index[i][0]
                task_name = Job_Data.index[i][1]
                k = mapper.get_job(job_name)
                required_data = Job_Data[data_partition_name][i]
                if not pd.isnull(required_data):
                    self.d_kis[k][i][s] = required_data

    def get_d_kis(self):
        return self.d_kis


class C_kij:
    def __init__(self, D_kis, mapper, Inter_Link):
        self.D_kis = D_kis
        self.mapper = mapper
        self.Inter_Link = Inter_Link
        self.datacenter_list = mapper.get_datacenter_list()

    def compute_c_kij(self):
        job_list = self.mapper.get_job_list()
        task_list = self.mapper.get_task_list()
        datacenter_list = self.mapper.get_datacenter_list()
        data_partition_list = self.mapper.get_data_partition_list()

        max_k = len(job_list)
        max_i = len(task_list)
        max_j = len(datacenter_list)
        max_s = len(datacenter_list)
        self.c_kij = np.zeros((max_k, max_i, max_j))
        d_kis = self.D_kis.get_d_kis()
        flag = np.zeros((max_k, max_i, max_j))  # 用来判断放在dc j是否无法获得所需的数据
        for k in range(max_k):
            for i in range(max_i):
                for s in range(max_s):
                    if d_kis[k][i][s] > 0:
                        for j in range(max_j):

                            # for s in range(max_s):
                            # if s != j and d_kis[k][i][s] > 0:  # `k`th task `i`th job need read data in datacenter `s`
                            d_kis_num = d_kis[k][i][s]
                            b_sj_num = self.b_sj[s][j]  # speed from `s` to `j`
                            if b_sj_num > 0 and flag[k][i][j] == 0:  # speed need > 0
                                cur_c = d_kis_num / b_sj_num
                                # max_c = max(max_c, cur_c)
                                # else:
                                #     max_c = 0  # there is no way from `s` to `j`, so we ignore `c[k][i][j]` and set it 0
                                #     break
                                self.c_kij[k][i][j] = max(self.c_kij[k][i][j], cur_c)
                            else:
                                flag[k][i][j] = 1
                                self.c_kij[k][i][j] = 0

    def compute_b_sj(self):  # checked: OK
        max_len = len(self.datacenter_list)
        self.b_sj = np.zeros((max_len, max_len))
        for s in range(len(self.datacenter_list)):
            for j in range(len(self.datacenter_list)):
                # if s != j:
                datacenter_s_name = self.mapper.get_datacenter(j)
                speed = self.Inter_Link[datacenter_s_name][s]
                if not pd.isnull(speed) and speed != '-':
                    self.b_sj[s][j] = speed

    def get_c_kij(self):
        return self.c_kij

    def get_b_sj(self):
        return self.b_sj


class E_kij:
    def __init__(self, Execution_Time, c_kij, Job_Task_List, mapper):
        self.e_kij = np.zeros(c_kij.shape)
        max_k = self.e_kij.shape[0]
        max_i = self.e_kij.shape[1]
        max_j = self.e_kij.shape[2]
        for n in range(len(Job_Task_List)):
            k_name = Job_Task_List[n][0]  # str job
            i_name = Job_Task_List[n][1]  # str task
            k = mapper.get_job(k_name)
            i = mapper.get_task(i_name)
            self.e_kij[k][i][0:max_j] = Execution_Time['Execution Time (s)'][i]

    def get_e_kij(self):
        return self.e_kij


class A_j:
    def __init__(self, Slot_Number):
        dc_num = len(Slot_Number['DC'])
        self.a_j = np.zeros(dc_num)
        for i in range(dc_num):
            self.a_j[i] = Slot_Number['Num of Slots'][i]

    def get_a_j(self):
        return self.a_j

    def decrease_a_j(self, index):
        self.a_j[index] = self.a_j[index] - 1


def get_M(mapper, Job_List, Job_Task_Dict):
    datacenter_list = mapper.get_datacenter_list()
    J = len(datacenter_list)
    sum_nk = 0
    for i in range(len(Job_List)):
        job_name = Job_List[i]
        sum_nk = sum_nk + len(Job_Task_Dict[job_name])
    return sum_nk * J


class Precedence:
    # i x i boolean matrix, the i th row task has several tasks of constraints
    def __init__(self, Task_List, mapper, Job_Precedence=None):
        self.Task_List = Task_List
        self.mapper = mapper
        if Job_Precedence is not None:
            self.Job_Precedence = Job_Precedence
        self.max_i = len(Task_List)
        self.precedence = np.zeros((self.max_i, self.max_i))

    def compute_precedence(self):
        str_constr_list = self.Job_Precedence['Precedence Constraint.1']
        for num in range(len(str_constr_list)):
            if len(str_constr_list[num]) > 3:  # { } == 3
                # the num th task has constraints
                for i in range(len(self.Task_List)):
                    if i != num and self.Task_List[i] in str_constr_list[num]:
                        self.precedence[num][i] = 1

    def get_precedence(self):
        self.compute_precedence()
        return self.precedence


class W_kij:
    def __init__(self, c_kij, e_kij, precedence):
        self.max_k, self.max_i, self.max_j = np.shape(c_kij)
        self.c_kij = c_kij
        self.e_kij = e_kij
        self.x_kij = np.ones((self.max_k, self.max_i, self.max_j))
        self.precedence = precedence

    def update_x_kij(self, x_kij):
        self.x_kij = x_kij

    def compute_w_kij(self):
        self.w_kij = np.zeros((self.max_k, self.max_i, self.max_j))
        for k in range(self.max_k):
            for i in range(self.max_i):
                max_w = 0
                for i_constr in range(self.max_i):
                    if self.precedence[i][i_constr] == 1:
                        for j in range(self.max_j):
                            if self.c_kij[k][i_constr][j] != 0:
                                cur_w = self.x_kij[k][i_constr][j] * (
                                        self.c_kij[k][i_constr][j] + self.e_kij[k][i_constr][j] +
                                        self.w_kij[k][i_constr][j])
                                if cur_w > max_w:
                                    max_w = cur_w
                for j in range(self.max_j):
                    self.w_kij[k][i][j] = max_w

    def get_w_kij(self):
        return self.w_kij

# if __name__ == "__main__":
#     # test code
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
#     Precedence_instance = Precedence(Task_List, mapper_instance, Job_Precedence=Job_Precedence)
#     precedence = Precedence_instance.get_precedence()
#
#     W_kij_instance = W_kij(c_kij=C_kij_instance.get_c_kij(), e_kij=E_kij_instance.get_e_kij(), precedence=precedence)
#     W_kij_instance.compute_w_kij()
#     w_kij = W_kij_instance.get_w_kij()
#     print("pause")
