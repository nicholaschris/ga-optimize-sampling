#!usr/bin/python
#ga.py
#genetic_algorithm_version5000.py

import datetime
import numpy as np
#from matplotlib import pyplot as plt
import make_population
import make_gene_map
import ga_class
from numpy import ma

from scipy.interpolate import Rbf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

my_gene_map = make_gene_map.GeneMap()
my_gene_map.get_array_attributes()
my_gene_map.make_gene_map_2()

bin_str_len = my_gene_map.string_length
population_array = my_gene_map.array

my_pop = make_population.InitPopulation()
my_pop.get_gene_size(bin_str_len)
my_pop.set_up_new_pop()

genetic = ga_class.GeneticAlgorithm()

genetic.population = my_pop.pop_list
genetic.gene_map = my_gene_map.gene_map
genetic.array = my_gene_map.array
genetic.lat_len = my_gene_map.lat_len
genetic.lon_len = my_gene_map.lon_len
genetic.string_length = my_gene_map.string_length


fittest_list = []
least_fit_list = []

gen_count = 0
while gen_count < 100:
	genetic.add_fitness_key()
	fittest_list = np.append(fittest_list, genetic.fittest_fitness)
	least_fit_list = np.append(least_fit_list, genetic.least_fitness)
	# print len(genetic.population)
	genetic.sort_pop()
	genetic.make_new_pop()
	genetic.elitism()
	genetic.population = genetic.new_pop
	gen_count+=1
	

dt = datetime.datetime.now()
date_time_str = dt.ctime()
date_time_short = date_time_str[11:]

print date_time_short

plt.close('all')
#plt.ylim(ymin=50, ymax=110)
plt.plot(fittest_list, label='Fittest Individuals')
plt.plot(least_fit_list, label='Least Fit Individuals')
plt.title('Genetic Algorithm - ' + str(date_time_short))
plt.xlabel('Generation')
plt.ylabel('fitness')
plt.legend(loc=7)
plt.savefig('/home/nicholas/figures/GA_test' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')

solution = genetic.population[0]['chrom_list']

def calc_chrom(chromosome, dictionary):
	gene_stepper = 0
	values_list = []
	coord_list = [(0, 0, 1)]
	gene_length = genetic.string_length -1 # hard code
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
	return coord_list, values_list

coord_list, values_list = calc_chrom(solution, genetic.gene_map)

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

# use RBF
rbf = Rbf(zs, ys, values_list, function='gaussian', epsilon=4)
ZI = rbf(loni, lati)

plt.close('all')
# plot the result
#n = plt.normalize(-2., 2.)
plt.subplot(1, 1, 1)
plt.pcolor(loni, lati, ZI, cmap=cm.jet, vmin=0, vmax=20)
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=0, vmax=20)
plt.title('RBF interpolation - multiquadrics')
#plt.xlim(0, 2)
#plt.ylim(0, 2)
plt.axis('tight')
plt.colorbar()
plt.savefig('rbf2d' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')

genetic.array = ma.masked_values(genetic.array, 1e20)
GI = genetic.array[0, :, :]
ZI = np.swapaxes(ZI, 1, 0)
diff_array = np.sqrt((GI - ZI)**2)
plt.pcolor(diff_array, cmap=cm.jet, vmin=0, vmax=10); 
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=0, vmax=20)
plt.axis('tight')
plt.colorbar()
plt.savefig('diff_array' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')
fitness = np.sqrt(np.sum((GI - ZI)**2))
print "The fitness of the solution is: %g" % fitness

plt.pcolor(GI, cmap=cm.jet, vmin=0, vmax=20)
plt.scatter(zs, ys, s, values_list, cmap=cm.jet, vmin=0, vmax=20)
plt.axis('tight')
plt.colorbar()
plt.savefig('snap_shot' + '-' + str(dt.year) + '-' + str(dt.month) + '-' + str(dt.day) + '-' + str(dt.hour) + 'h' + str(dt.minute) + '.png')
plt.close('all')