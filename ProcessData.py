import numpy as np
from utils import TwoWayDict
import pandas as pd
from algorithm_project import get_data, job_task_relationship, get_index
import math

def get_K(Job_List):
    return len(Job_List)

def get_D():
    pass

def get_T_k(Job_Task_List, K):
    T_k = {}
    for i in range(K):
        T_k[K] = len(Job_Task_List['A'+i])

class Mapper:
    def __init__(self, Job_List, Task_List, Slot_Number, Data_Partition):
        self.job = TwoWayDict()
        self.task = TwoWayDict()
        self.datacenter = TwoWayDict()
        self.slots = dict() # 每个datacenter有几个slot
        self.data_partition2datacenter = dict()
        #self.datacenter2data_partition = dict()
        for i in range(len(Job_List)):
            self.job[Job_List[i]]=i
        for i in range(len(Task_List)):
            self.task[Task_List[i]]=i
        for i in range(len(Slot_Number['DC'])):
            self.datacenter[Slot_Number['DC'][i]]=i
        for i in range(len(Slot_Number['Num of Slots'])):
            self.slots[i] = Slot_Number['Num of Slots'][i] # 例，0对应2
            self.slots[self.datacenter[i]] = Slot_Number['Num of Slots'][i] # 例，DC1对应2
        for i in range(len(Data_Partition['Data Partition'])):
            data_partition_name = Data_Partition['Data Partition'][i]
            datacenter_name = Data_Partition['Location'][i]
            self.data_partition2datacenter[data_partition_name] = datacenter_name
            #self.datacenter2data_partition[datacenter_name] = data_partition_name

        self.job_list = Job_List
        self.task_list = Task_List
        self.datacenter_list = list()
        for i in range(len(Slot_Number['DC'])):
            self.datacenter_list.append(Slot_Number['DC'][i])
        self.data_partition_list = list()
        for i in range(len(Data_Partition['Data Partition'])):
            self.data_partition_list.append(Data_Partition['Data Partition'][i])


    def get_job(self,key):
        return self.job[key]
    def get_task(self,key):
        return self.task[key]
    def get_datacenter(self,key):
        return self.datacenter[key]
    def get_slots(self,key):
        return self.slots[key]
    def get_data_partition2data_center(self,key):
        return self.data_partition2datacenter[key]

    def get_job_list(self):
        return self.job_list
    def get_task_list(self):
        return self.task_list
    def get_datacenter_list(self):
        return self.datacenter_list
    def get_data_partition_list(self):
        return self.data_partition_list



class D_kis: # checked: OK
    def __init__(self, mapper, Job_Data):
        self.mapper = mapper
        job_list = mapper.get_job_list()
        task_list = mapper.get_task_list()
        datacenter_list = mapper.get_datacenter_list()
        data_partition_list = mapper.get_data_partition_list()
        max_k = len(job_list)
        max_i = len(task_list)
        #max_s = len(data_partition_list)
        max_s = len(datacenter_list)

        max_datapartition = len(data_partition_list)


        # assign task_id to d[k][i][s]
        self.d_kis = np.zeros((max_k,max_i,max_s))
        #for k in range(max_k):
            #for i in range(max_i):
        for ss in range(max_datapartition):
            data_partition_name = data_partition_list[ss]
            datacenter_name = mapper.get_data_partition2data_center(data_partition_name)
            s = mapper.get_datacenter(datacenter_name)
            for i in range(max_i):
                job_name = Job_Data.index[i][0]
                task_name = Job_Data.index[i][1]
                k=mapper.get_job(job_name)
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
        self.c_kij = np.zeros((max_k,max_i,max_j))
        d_kis = self.D_kis.get_d_kis()
        for k in range(max_k):
            for i in range(max_i):
                for j in range(max_j):
                    max_c = 0
                    for s in range(max_s):
                        if s!=j and d_kis[k][i][s] > 0: # `k`th task `i`th job need read data in datacenter `s`
                            d_kis_num = d_kis[k][i][s]
                            b_sj_num = self.b_sj[s][j] # speed from `s` to `j`
                            if b_sj_num > 0: # speed need > 0
                                cur_c = d_kis_num / b_sj_num
                                max_c = max(max_c, cur_c)
                            else:
                                max_c = 0 # there is no way from `s` to `j`, so we ignore `c[k][i][j]` and set it 0
                                break
                    self.c_kij[k][i][j] = max_c
        # for i in range(max_i):
        #     for j in range(max_j):
        #         for k in range(max_k):
        #             max_c = -1
        #             for s in range(max_s):
        #                 if s != j:
        #                     d_kis_num = d_kis[k][i][s]
        #                     b_sj_num = self.b_sj[s][j]
        #                     if b_sj_num != 0:
        #                         cur_c = d_kis_num / b_sj_num
        #                         max_c = max(max_c, cur_c)
        #             self.c_kij[k][i][j] = max_c

    def compute_b_sj(self): # checked: OK
        max_len = len(self.datacenter_list)
        self.b_sj = np.zeros((max_len,max_len))
        for s in range(len(self.datacenter_list)):
            for j in range(len(self.datacenter_list)):
                if s!=j:
                    datacenter_s_name = self.mapper.get_datacenter(j)
                    speed = self.Inter_Link[datacenter_s_name][s]
                    if not pd.isnull(speed) and speed != '-':
                        self.b_sj[s][j] = speed

class c_kij_e_kij:
    def __init__(self, Job_List, Job_Task_Dict):
        pass


if __name__ =="__main__":
    # test code
    Inter_Link, Job_Data, Execution_Time, Job_Precedence, Slot_Number, Data_Partition = get_data()
    Job_List, Task_List, Job_Task_List, Job_Task_Dict = job_task_relationship()
    get_index(Job_Data, Job_Task_List)

    mapper_instance = Mapper(Job_List, Task_List, Slot_Number, Data_Partition)
    D_kis_instance = D_kis(mapper_instance, Job_Data)
    C_kij_instance = C_kij(D_kis_instance,mapper_instance,Inter_Link)
    C_kij_instance.compute_b_sj()
    C_kij_instance.compute_c_kij()

    print("pause")
    print("pause")



