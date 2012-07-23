# calc_fitness.py

''' This module calculates various values (mean, stdev) of a chromosome's value list.
To add: after calculating fitness a fitness value shoud be assigned to chromosome.

'''

import numpy as np
from numpy import ma

def calc_fitness(values_list, pop_array):
	print 'The mean value of the sample is %g.' % ma.mean(ma.masked_values(values_list, 1e+20))
	print 'The stdev of the sample is %g.' % ma.std(ma.masked_values(values_list, 1e20))
	print "The count of non_masked values is %g." % ma.count(ma.masked_values(values_list, 1e20))
	# fitness function
	fitness = np.abs(ma.mean(ma.masked_values(values_list, 1e+20)) - ma.mean(pop_array))\
	               + np.abs(ma.std(ma.masked_values(values_list, 1e20)) - ma.std(pop_array))
	fitness = 200 - fitness*ma.count_masked(ma.masked_values(values_list, 1e20))
	print "The fitness of the chromosome/solution is %g." % fitness
	return fitness
