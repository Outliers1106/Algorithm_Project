# -*- coding: utf-8 -*-
"""ModifySimulationData.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1mGU3Ds8UjPNCbdI4uuq6VkjBC23VOYG4
"""

import pandas as pd
import numpy as np
import random

def read_data():
  Inter_Link=pd.read_excel('SimulationData.xlsx',sheet_name='Inter-DCL')
  Execution_Time=pd.read_excel('SimulationData.xlsx',sheet_name='JobList',usecols=[845])
  Slot_Number=pd.read_excel('SimulationData.xlsx',sheet_name='DC-Details',usecols=[0,1])
  Job_Data=pd.read_excel('SimulationData.xlsx',sheet_name='JobList',usecols=range(2,845))
  return Inter_Link, Execution_Time, Slot_Number, Job_Data

def InterLinkModify(Inter_Link):
  for i in range(len(Inter_Link)):
    Inter_Link.iloc[i,i+1]=1000
  for i in range(227):
    for j in range(1,228):
      if Inter_Link.iloc[i,j]<10:
        Inter_Link.iloc[i,j]=10
  return Inter_Link

def JobDataModify(Job_Data):
  task=pd.read_excel('SimulationData.xlsx',sheet_name='JobList',usecols=[1])
  task_list=[]                                               
  for i in range(len(task)):
    task_list.append(task['Unnamed: 1'].iloc[i])  
  Job_Data.index=task_list
  Precedence=pd.DataFrame(index=task_list,columns=['Precedence Constraint'])
  Precedence=Precedence.fillna('{}')
  for i in range(566): 
    for j in range(843):  #遍历表格中的每一个值
      if not Job_Data.index[i][1] in Job_Data.columns[j]:  #每一个task不会从其他的job中获取数据
        Job_Data.iloc[i,j]=0
      else:
        if not Job_Data.columns[j][0] in 't':  #有些task可能会从Data Partition（例如A1,A2等）获取数据
          p=random.uniform(0,1)
          if p>=0.5:                    #产生需要的概率为0.5
            Job_Data.iloc[i,j]=0
        elif Job_Data.columns[j][0] in 't' and int(Job_Data.index[i][2:])>int(Job_Data.columns[j][2:]):   #可能产生Precedence的组合                      
          q=random.uniform(0,1)
          if q>=0.4:                    #产生precedence constraint的概率为0.4
            Job_Data.iloc[i,j]=0
          else:
            Precedence.iloc[i,0]=Precedence.iloc[i,0][:-1]+str((Job_Data.columns[j],Job_Data.index[i]))+',}'
        else:                           #不可能产生precedence constraint的task组合，也不会有数据传输
          Job_Data.iloc[i,j]=0
  return Job_Data, Precedence

def DataPartitionFinish(Job_Data):
  Partition_list=[]
  for j in range(843):
    if not Job_Data.columns[j][0] in 't':
      Partition_list.append(Job_Data.columns[j])
  Data_Partition=pd.DataFrame(index=Partition_list,columns=['Location'])
  for i in range(len(Data_Partition)):
    location='DC'+str(random.randint(1,227))
    Data_Partition.iloc[i][0]=location
  return Data_Partition

def SaveData(Inter_Link,Job_Data, Execution_Time, Precedence, Slot_Number, Data_Partition):
  writer = pd.ExcelWriter('SimulationPython.xlsx')
  Inter_Link.to_excel(writer,sheet_name='Inter-Datacenter Links',index=False)
  Job_Data.to_excel(writer,sheet_name='Job List')
  Execution_Time.to_excel(writer,sheet_name='Execution Time')
  Precedence.to_excel(writer,sheet_name='Precedence')
  Slot_Number.to_excel(writer,sheet_name='Slot Number',index=False)
  Data_Partition.to_excel(writer,sheet_name='Data Partition')
  writer.save()
  writer.close()

def main():
  Inter_Link, Execution_Time, Slot_Number, Job_Data=read_data()
  Inter_Link=InterLinkModify(Inter_Link)
  Job_Data,Precedence=JobDataModify(Job_Data)
  Data_Partition=DataPartitionFinish(Job_Data)
  SaveData(Inter_Link,Job_Data, Execution_Time, Precedence, Slot_Number, Data_Partition)

main()

