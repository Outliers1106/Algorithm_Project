from ReadData import get_data, job_task_relationship, get_index
import gurobipy

if __name__ == "__main__":
    Inter_Link, Job_Data, Execution_Time, Job_Precedence, Slot_Number, Data_Partition = get_data()
    Job_List,Task_List,Job_Task_List, Job_Task_Dict = job_task_relationship()
    get_index(Job_Data,Job_Task_List)
    print("a")