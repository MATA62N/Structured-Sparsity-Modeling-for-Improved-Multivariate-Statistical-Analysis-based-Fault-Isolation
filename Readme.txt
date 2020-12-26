#Rawdata.mat is a simulation data set obtained according to the data generation method mentioned in Section 5. The code realizes all the results of the multiplicative fault for the simulation experiments in Section 5, and the diagnosis methods of sensor bias and coal pulverizer are similar.

# Usage
1.Run the script *runsim.m*
```
2.Subroutine: Fault isolation based on Group Lasso
```	
./simulate_GL.m
```	
3.Subroutine: Fault isolation based on Sparse Group Lasso
```
./simulate_SGL.m
```
4.Subroutine: Fault isolation based on tree-structured sparsity
```
./simulate_tree.m
```
5.Subroutine: Fault isolation based on partially known sparse support
```
./simulate_know.m
```
6.Subroutine: Fault isolation based on conventional reconstruction
```
./simulate_simple.m
```
7.Subroutine: Fault isolation based on l_1 penalty
```
./simulate_l1.m
```