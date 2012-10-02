#!usr/bin/python
#ga.py
#genetic_algorithm_version5000.py
# 	TODO
#	add comments
#	make sure it works
#	check why the algorithm never chooses points above lat=35
#	why does the final fitness message differ from the loop fitness message?
#	What is 4095 and what is 7200?
#	Does sort pop work if it is a list and not dictionary?
#	Using inheritance or other OOP thingies?

import datetime
import numpy as np
import make_population
import make_gene_map
import matplotlib
import ga_class
import matplotlib.pyplot as plt
from numpy import ma
from scipy.interpolate import Rbf
from matplotlib import cm

matplotlib.use('Agg')

my_gene_map = make_gene_map.GeneMap()	# creates a new instrance of the GeneMap class
my_gene_map.get_array_attributes()			# gets the attributes of the array that was used to create the gene_map
my_gene_map.make_gene_map_2()			# creates a gene_map using the 2nd make_gene_map method

print np.shape(my_gene_map.array) 
print my_gene_map.array_size
print my_gene_map.string_length

bin_str_len = my_gene_map.string_length		# get the length of the bit string
population_array = my_gene_map.array               # Where is this used?

my_pop = make_population.InitPopulation()	# creates a new instrance of the InitPopulation class
my_pop.get_gene_size(bin_str_len)			# gets the size of the gene from the string length of genes in genemap
my_pop.set_up_new_pop()					# sets up a new population from the 

genetic = ga_class.GeneticAlgorithm()			# creates an instance of the genetic algorithm class

genetic.population = my_pop.pop_list			# the population for the genetic algorithm is the same as the my_pop. A list of dictionaries.
genetic.gene_map = my_gene_map.gene_map	# the genemapfor the GA taken from genemap class
genetic.array = my_gene_map.array
genetic.time_len = my_gene_map.time_len#  array from genemap
genetic.lat_len = my_gene_map.lat_len		# other attributes from genemap
genetic.lon_len = my_gene_map.lon_len
genetic.string_length = my_gene_map.string_length


fittest_list = []	# Use the fittest list to plot evolution of algorithm DOES IT WORK?
least_fit_list = []	# likewise for least fit list

NUM_GEN = 1000

gen_count = 0
while gen_count < NUM_GEN:
	genetic.add_fitness_key()									# calculate fitness for each chromosome in population and add fitness key
	fittest_list = np.append(fittest_list, genetic.fittest_fitness)		# use to plot progress of algorithm
	least_fit_list = np.append(least_fit_list, genetic.least_fitness)	# use to plot progress of algorithm
	#~ print len(genetic.population)
	genetic.sort_pop()											# sorts the population from least fit to fittest - DOES IT WORK?
	genetic.make_new_pop()										# Uses selection, crossover, mutation to make a new population
	genetic.elitism()												# Replaces a random chromosome with the fittest chromosome from last pop SHOULD IT BE USED
	genetic.population = genetic.new_pop
	# Replaces the old population with the new popualtion so it can be run again
	print "There are %g generations left" %(NUM_GEN-gen_count) 
	gen_count+=1												# the counter
	
#~ genetic.add_fitness_key()

### Some stuff to save time and date on plot and filenames
dt = datetime.datetime.now()
date_time_str = dt.ctime()
date_time_short = date_time_str[11:]
print date_time_short

plt.close('all')

### A plot to visualise how the GA evolved ###
plt.plot(fittest_list, label='Fittest Individuals')
plt.plot(least_fit_list, label='Least Fit Individuals')
plt.title('Genetic Algorithm - ' + str(date_time_short))
plt.xlabel('Generation')
plt.ylabel('fitness')
plt.legend(loc=7)
plt.savefig('/home/nicholas/figures/GA_test' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
### Ebd plotting ###

# Does the population not need to be sorted first? Huh???? 
genetic.add_fitness_key()
genetic.sort_pop()
# Does the population not need to be sorted first? Huh???? 
solution = genetic.population[0]['chrom_list']					# This should be the solution (the best locations to sample at)
# Does the population not need to be sorted first? Huh???? 

### This saves the coordinate tuples in a list ###
# But do I need the values too? Yes for the RBF...
# The list should hold 200 coordinates otherwise something is wrong.
def calc_chrom(chromosome, dictionary):
	'''
	This saves the coordinate tuples in a list.
	'''
	gene_stepper = 0
	values_list = []
	coord_list = [(0, 0, 1)]
	gene_length = genetic.string_length  # hard code REMOVEDD MINUS ONE
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
	print np.size(coord_list) # DEBUG
	return coord_list, values_list

coord_list, values_list = calc_chrom(solution, genetic.gene_map)

# WHY DOES COORD_LIST HAVE 216 COORDINATES AND NOT 200 [solved?]

#values_list = ma.masked_where(values_list=0)

coords_vals = dict(zip(coord_list, values_list))
if coords_vals.has_key((999, 999, 999)):	
	coords_vals.pop((999,999,999))
	
x_nodes = []	
y_nodes = []
z_nodes = []
values_list = []
for item in coords_vals:
	x_nodes.append(item[0])
	y_nodes.append(item[1])
	z_nodes.append(item[2])
	values_list.append(coords_vals[item])
	
	
'''coord_list = np.unique(coords_vals)
values_list = coords_vals.values()'''

xs = x_nodes
ys = y_nodes
zs = z_nodes
values_list = np.array(values_list)

lat = my_gene_map.lat_len
lon = my_gene_map.lon_len

ti_lat = np.arange(0, lat)
ti_lon = np.arange(0, lon)

lati, loni = np.meshgrid(ti_lat, ti_lon)

s = len(xs)


from Scientific.IO.NetCDF import NetCDFFile
### get landmask
nc = NetCDFFile('/home/nicholas/data/ORCA2_landmask.nc','r')
mask = ma.masked_values(nc.variables['MASK'][:, :, :40, :180], -9.99999979e+33)
nc.close()

### Trying out 3d RBF (insert into ga_class)
### meshgrid 3d
xxx, yyy, zzz = np.lib.index_tricks.mgrid[0:10, 0:40, 0:180]


# use RBF
rbf = Rbf(xs, ys, zs, values_list, function='gaussian', epsilon=4)
#~ rbf = Rbf(zs, ys, values_list, function='gaussian', epsilon=4)
ZI = rbf(xxx, yyy, zzz)
#~ ZI = rbf(loni, lati)
#~ ZI = np.swapaxes(ZI, 1, 0)
ZI = ZI**mask[0, 0:10, :, :]

plt.close('all')

plt.subplot(1, 1, 1)
plt.pcolormesh(ma.mean(ZI, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.title('RBF interpolation - multiquadrics')

plt.axis('tight')
plt.colorbar()
plt.savefig('rbf2d' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')

genetic.array = ma.masked_values(genetic.array, 1e20)
GI = genetic.array[:, :, :]
print np.shape(GI) # DEBUG

diff_array = np.sqrt((GI - ZI)**2)

plt.pcolormesh(ma.mean(diff_array, axis=0), cmap=cm.jet, vmin=0, vmax=10); 
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.axis('tight')
plt.colorbar()
plt.savefig('diff_array' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')
fitness = np.sqrt(np.sum((GI - ZI)**2))
print np.mean(GI), np.mean(ZI), np.std(GI), np.std(ZI)
print "The fitness of the solution is: %g" % fitness # This answer should be the same as the final fittest solution???

plt.pcolormesh(ma.mean(GI, axis=0), cmap=cm.jet, vmin=-5, vmax=10)
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=-5, vmax=10)
plt.axis('tight')
plt.colorbar()
plt.savefig('snap_shot' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')



