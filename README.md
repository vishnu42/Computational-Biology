# CompBio
Computational Biology Final Project:

# Generating graphs 
New Signature scheme for alignment free network comparison of
>  python data_generator_rewire.py --N 1000 --rho 0.01 --num 3 --seed 10 --rewire 10

will create a folder tree containing 10 graphs for each of the four models with
a edge density of 0.1% and 1000 nodes. It also generates 1 folder with 10% rewired
edges.


# Running covariance_matrix.py: (basic version)
1) The program expects an input graphs dataset, and output file into which the result will be written and a parameter k(number of power iterations, default = 5).
2) command to run the program:	python covariance_matrix.py --input dataset_3 --output outfile --k 5

# Running covariance_matrix_threads.py: (Multi threading, dynamic programming)
1) The program expects an input graphs dataset, and output file into which the result will be written, a parameter k(number of power iterations, default = 5) and no_of_threads (default 10).
2) command to run the program:	python covariance_matrix_threads.py --input dataset_3 --output outfile --k 5 --threads 5
