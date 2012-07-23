import datetime
import numpy as np
#from matplotlib import pyplot as plt
import make_population
import make_gene_map
import ga_class
from numpy import ma

from scipy.interpolate import Rbf
'''import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm'''

def calc_chrom_2(self, index):
	dictionary = self.gene_map
	chromosome = self.population[index]['chrom_list']
	gene_stepper = 0
	values_list = []
	coord_list = [(0, 0, 1)]
	gene_length = 11 # hard code
	while gene_stepper < len(chromosome) - gene_length:
		gene_list = chromosome[gene_stepper:gene_stepper+20]
		current_gene = ''
		bit_stepper = 0
		while bit_stepper < gene_length:
			current_gene = current_gene + gene_list[bit_stepper]
			bit_stepper += 1
		# print 'Current gene is %s.' % current_gene
		values_list = np.append(values_list, dictionary[current_gene]['value'])
		coord_tuple = dictionary[current_gene]['coordinate']
		# print coord_tuple # debug
		# print coord_tuple # debug
		coord_list.append(coord_tuple)
		gene_stepper += gene_length
	coord_list.pop(0)
	#return coord_list, values_list

	coord_list, values_list = calc_chrom(solution, genetic.gene_map)

	coords_vals = dict(zip(coord_list, values_list))
	coords_vals.pop((999,999,999))
	coord_list = np.unique(coords_vals)
	values_list = coords_vals.values()

	# print coord_list


	xs = coord_list[:, 0]
	ys = coord_list[:, 1]
	zs = coord_list[:, 2]
	values_list = np.array(values_list)

	lat = my_gene_map.lat_len
	lon = my_gene_map.lon_len

	ti_lat = np.arange(0, lat)
	ti_lon = np.arange(0, lon)

	lati, loni = np.meshgrid(ti_lat, ti_lon)

	s = len(xs)

	# use RBF
	rbf = Rbf(zs, ys, values_list, epsilon=2)
	ZI = rbf(loni, lati)

	GI = genetic.array[0, :, :]
	GI = np.swapaxes(GI, 1, 0)
	diff_array = (GI - ZI)**2
	fitness = np.sum((GI - ZI)**2)
	print "The fitness of the solution is: %g" % fitness
	self.population[index]['fitness'] = fitness