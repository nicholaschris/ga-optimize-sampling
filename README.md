ga-optimize-sampling

A Genetic Algorithm implementation in python that optimizes
sampling in 3D space.

In load\_cflux\_masked.py:
	A variable is selected from a netcdf file which is saved as
	a numpy array. A mask can be applied and the indices chosen.

In make\_gene\_map.py:
	The numpy array is used to create a gene map for the 
	genetic algorithm.
	Each data point in the array is assigned a binary string.
	The binary string is known as the gene string.
	A dictionary is created where each binary string is the key.
	The values of these keys are also dictionaries with the 
	keys "location" and "value" and their associated values.

	eg.{'001101': {'location': (0, 10, 5), 'value': 1.2},
	    '001111': {'location': (0, 19, 6), 'value': 3.0}} 

	If the last binary string is has zeros in it eg '101111'
	then more binary strings are added to the dictionary
	until all the possible permuations of 0, 1 are used.
	These value and loaction for these binary strings are
	1e09 and (999, 999, 999). 

In make\_population.py:
	The initial population is created. 
	This is done by randomly selecting a zero or a one and 
	appending it to a chromosome string.
	The length of the chromosome string is the product of the length
	of a gene string and the number of locations you would like 
	to select to use as a sample.
	This process is then repeated N times with N being the size
	you would like the population to be.

In ga\_class.py:
	This is where you will find the attributes and methods that you
	can use for the genetic algorithm such as selection, crossover, 
	mutation, and the fitness function.

In ga.py:
	This file runs the genetic algorithms N times. Once the 
	optimal sampling solution has been determined by the ga
	the results will be plotted.
	 
