[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An efficient branch-and-price algorithm for the submodular bin packing problem

This repository accompanies the paper submitted to the [INFORMS Journal on Computing (IJOC)](https://pubsonline.informs.org/journal/ijoc).
 All source code is released under the [MIT License](LICENSE).

## 1. Citation

If you use this repository or the accompanying paper in your work, please cite both the paper and this repository using their respective DOIs.

<!-- TODO: Add DOI of the paper -->
<!-- TODO: Add DOI of this repository -->

BibTeX for citing this snapshot of the repository:

<!-- TODO: Insert BibTeX here -->

## 2. Description

This repository provides the implementation of our algorithms to reproduce the computational results reported in the paper. The software allows users to solve both the **Submodular Knapsack Problem (SMKP)** and the **Submodular Bin Packing Problem (SMBP)**

## 2.1 Building

The code is designed for **Linux** environments. 

1. Clone the repository:

```bash
git clone <repo-link>
cd <repo-folder>
```

2. Build the code using **CMake**:

```bash
mkdir build && cd build
cmake ..
make
```

This will generate the executable `SMKMP` in the `build` directory. Since we use Gurobi as the commercial solver, the users need to change the path related to the Gurobi solver in the CMakeLists.txt.

## 2.2 Running the Software

 The executable requires command-line arguments. There are **8 parameters in total**, three mandatory and five optional:

```makefile
-i   <file>   # input .txt file (mandatory)
-o   <file>   # output .txt file (mandatory)
-p   <SMKP | SMBP>  # problem type (mandatory)
-m   <Gurobi>   # solution method (optional, default: DP for SMKP, BPC for SMBP)
-a   <ULA-VD | ULA-FD | LLA-VD | LLA-FD>  # algorithm combinations for SMBP (optional, default: ULA-VD)
-S            # toggle SR3 inequalities (optional, default: use SR3; include flag to disable)
-D            # toggle DIs (optional, default: not used; include flag to enable)
-r            # required for instances from Xu et al. (2023); ignored otherwise

```

The "-i", "-o", and "-p" are mandatory parameters, "-m", "-a", "-S", and "-D" are optional parameters, "-r" is necessary for instances from Xu et al.(2023) and useless for other instances.

## 2.3 Example Commands

Solve SMBP (Xu et al., 2023 instances):

```bash
cd build
./SMKMP -i path/to/input.txt -o path/to/output.txt -p SMBP -r
```

Solve SMKP:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMKP
```

Solve SMBP without SR3 and with DIs:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -S -D -r
```

Call the Gurobi solver for SMKP:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMKP -m Gurobi
```

Call the Gurobi solver for SMBP:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -m Gurobi -r
```

Specify different algorithm combinations for SMBP:

```bash
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -a ULA-VD -r
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -a ULA-FD -r
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -a LLA-VD -r
cd build
./SMKMP -i path/to/the/input/file -o path/to/the/output/file -p SMBP -a LLA-FD -r
```



## 2.4 Example outputs

**SMKP with DP**

```
Solve SMKP use DP!

heuristic labeling get lower bound: 53
heuristic labeling use time: 0s
The best solution obtained by my labelsetting algorithm is: 53
The solution obtained by TS heuristic is: 53
The solution obtained by heuristic labeling algorithm is: 53
The time for my labelsetting algorithm is: 0.008s
The time for heuristic primal solution is: 0.003s
The time for linear relaxation (dual bound) is: 0s
The time for heuristic DP is: 0s
The time for exact dual bound is: 0s
The time for final exact DP is: 0s
The total number of generated labels is: 1948
The number of non-dominated labels is: 5
The number of CB fathomed labels is: 233
The number of dominated labels is: 740
The best item set is: 0,1,2,3,5,6,9,10,12,13,16,20,22,23,25,26,27,28,31,34,35,37,40,44,45,46,47,48,49,51,52,56,59,60,62,63,66,67,68,71,73,74,75,78,81,82,83,86,87,88,89,96,98
```
**SMKP with Gurobi**

```
Solve SMKP use Gurobi!
Optimal solution found (tolerance 1.00e-06)
Best objective 5.300000000000e+01, best bound 5.300000000000e+01, gap 0.0000%
***************The results obtained by Gurobi**************
The optimal objective value is :53
The solution time is :0.195881
The chosen items are as follows:
0,1,2,3,5,6,9,10,12,13,16,20,22,25,26,27,28,31,34,35,37,40,44,45,46,47,48,49,51,52,56,59,60,62,63,66,67,68,71,73,74,75,78,81,82,83,86,87,88,89,92,96,98,
The total number of chosen items is: 53
t = 8.02946
```
**SMBP with BPC (default: ULA-VD + SR3, without DIs)**

```
Solve SMBP use BPC!
Use algorithm combination: ULA-VD

The total number of generated labels is: 4169646
The number of non-dominated labels is: 250
The number of CB fathomed labels is: 28427
The number of dominated labels is: 2049238
Premature encountered, CG is stopped!
Heuristic DP time is: 1.62 s
Exact DP time is: 33.26 s
CG time is: 36.565 s
The root node CG iteration number is: 181
The number of DOIs added to the root nodeis: 0
The DOIs are NOT selected in the basis!
No branches is needed and the best solution use 6 bins!
The dual bound of the root node is 5.80641
The worst dual bound during branch 5.80641
The upper bound obtain by the HGS is 6 bins!
The best upper bound is 6 bins!
The total time for BP is 39.136
The total CG iteration number is: 181
The time spent on heuristic is 2.551
The time for solving root node is 36.566
The time spent on adjusting stages is 0
The time spent on pricing is 34.88
The average time for each node is 39.136
The number of pattern number branches are 0
The number of arc branches are 0
The number of nodes are 1
The maximum depth of the B_B tree is 0
The average depth of the B_B tree is 0
The number of columns generated by the HLA is 14355
The number of columns generated by the ELA is 58
```
**SMBP with Gurobi**

```
Solve SMBP use Gurobi!

***************The results obtained by Gurobi**************
The computational time is: 1.75089
The number of bins used is (upper bound): 6
Best Bound: 6
MIP Gap: 0
```

