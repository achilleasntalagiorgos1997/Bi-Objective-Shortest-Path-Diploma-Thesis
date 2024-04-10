# Introduction
This repository contains the code and documentation for my diploma thesis conducted at the University of Ioannina. The thesis focuses on the Biobjective Shortest Path Problem and presents a comparison between the Biobjective A* and Path-Pair-A* algorithms. The algorithms are implemented in C, a low-level programming language, to ensure efficiency and effectiveness.

# Contents
source/: This directory contains the source code of the implemented algorithms.
datasets/: Example datasets used for testing and evaluation purposes are stored here.
thesis/: The thesis document in Greek language.
Power-Point-Presentation/: The PP presentation in greek. (Both PDF and PP files).

# Biobjective Shortest Path Problem
The Biobjective Shortest Path Problem (BSP) involves finding the Pareto optimal set of solutions for two conflicting objectives: minimizing cost and minimizing distance. This problem has various real-world applications in transportation, logistics, and network optimization.

# Implemented Algorithms
1. Biobjective A* (BOA*)
BOA* is an adaptation of the classical A* algorithm to handle biobjective optimization problems. It explores the search space efficiently while maintaining a set of non-dominated solutions.
2. Path-Pair-A* (PPA*)
PPA* is another variant of the A* algorithm designed specifically for biobjective shortest path problems. It utilizes the concept of path pairs to efficiently compute the Pareto optimal set.


Diploma Thesis Repository
Introduction
This repository contains the code and documentation for my diploma thesis conducted at the University of Ioannina. The thesis focuses on the Biobjective Shortest Path Problem and presents a comparison between the Biobjective A* and Path-Pair-A* algorithms. The algorithms are implemented in C, a low-level programming language, to ensure efficiency and effectiveness.

Contents
source/: This directory contains the source code of the implemented algorithms.
datasets/: Example datasets used for testing and evaluation purposes are stored here.
thesis.pdf: The thesis document in Greek language.
Biobjective Shortest Path Problem
The Biobjective Shortest Path Problem (BSP) involves finding the Pareto optimal set of solutions for two conflicting objectives: minimizing cost and minimizing distance. This problem has various real-world applications in transportation, logistics, and network optimization.

Implemented Algorithms
1. Biobjective A* (Bi-A*)
Bi-A* is an adaptation of the classical A* algorithm to handle biobjective optimization problems. It explores the search space efficiently while maintaining a set of non-dominated solutions.

2. Path-Pair-A* (PPA*)
PPA* is another variant of the A* algorithm designed specifically for biobjective shortest path problems. It utilizes the concept of path pairs to efficiently compute the Pareto optimal set.

# Usage
To use the implemented algorithms, follow these steps:

1. Clone the repository to your local machine.

2. Navigate to the source/ directory.

3. Compile the source code using your preferred C compiler :
gcc -o biobjective_astar BOA.c
gcc -o path_pair_astar PPA.c

4. Execute the compiled programs with appropriate dataset inputs:
./biobjective_astar <dataset_file>
./path_pair_astar <dataset_file>

# Contributing
Contributions to this project are welcome. If you find any issues or have suggestions for improvements, please feel free to open an issue or create a pull request.


