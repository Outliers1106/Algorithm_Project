# Algorithm_Project
This is the group project for class Algorithm Design and Analysis (CS7310) in SJTU.

Team mates: Tu Yanlun , Xie Yuzhang, Yu Jingwei (Team 16)

## Requirements:
python==3.8.0
pandas==1.1.0 \
xlrd==1.2.0 \
numpy==1.19.1 \
gurobipy==9.1.0

You can install gurobipy by
```
conda config --add channels "http://conda.anaconda.org/gurobi"
conda install gurobi 
```
Make sure that you can access license of gurobi. 

## Usage
There are two folds, one for toy data, the other for simulated data.
```python
python SimData/main.py 
python ToyData/main.py 
``` 
After running, you can get the final assignment of matrix $x$ with numpy format.