#! usr/bin/python
# ga_class.py

'''Doc String.'''

import numpy as np
import random
from numpy import ma
from operator import itemgetter
import new_calc_chrom_for_rbf
from scipy.interpolate import Rbf

class GeneticAlgorithm(object):

	def __init__(self):
		self.population = []    	#get from Population class # change to list
		self.gene_map = {}      	# get from GeneMap class (gene_length too??)
		self.sum_fitness = 0
		self.new_pop = []
		self.array = []
		self.least_fitness = 0
		self.fittest_fitness = 0
		self.mutation_rate = 0.0000001
		self.diff_array = []
		self.string_length = 0
		self.fitness_mean = 100
		self.fitness_stdev = 100
		self.fitness_masked_count = 200
		
	def calc_chrom_rbf_part_1(self, index):
		
		'''
		Calculates the fitness of a chromosome according to 
		the square of the difference of an RBF interpolation of
		the sampled points from the chromosome and the
		model data.
		Used in add_fitness_key method.
		'''
		dictionary = self.gene_map
		chromosome = self.population[index]['chrom_list']
		gene_stepper = 0
		values_list = []
		coord_list = []
		gene_length = self.string_length # hard code ? AGAIN -1 ???
		while gene_stepper < len(chromosome) - gene_length: # The gene stepper goes from gene to gene o chromosome
			gene_list = chromosome[gene_stepper:gene_stepper+20] # chooses a gene
			current_gene = '' # initialises current gene
			bit_stepper = 0 # counter for the bits
			while bit_stepper < gene_length:
				current_gene = current_gene + gene_list[bit_stepper] #constructs the current gene
				bit_stepper += 1
			values_list = np.append(values_list, dictionary[current_gene]['value']) #appends the value of the current gene to a list
			coord_tuple = dictionary[current_gene]['coordinate'] # assigns  current coordinate to a variable
			#~ print coord_tuple
			coord_list.append(coord_tuple)
			#~ coord_list = np.append(coord_list, coord_tuple) # appends the coordinate of the current gene to the list
			#~ print coord_list 
			gene_stepper += gene_length # goes to the next gene in the loop
		
		coords_vals = dict(zip(coord_list, values_list))
		return coords_vals
		
	def 	calc_chrom_rbf_part_2(self, index, coords_vals):
		
		from Scientific.IO.NetCDF import NetCDFFile
		
		### get landmask
		nc = NetCDFFile('/home/nicholas/data/ORCA2_landmask.nc','r')
		mask = ma.masked_values(nc.variables['MASK'][:, :, :40, :180], -9.99999979e+33)
		nc.close()
		
		#~ if coords_vals.has_key((999, 999, 999)):	
			#~ coords_vals.pop((999,999,999))
		x_nodes = []	
		y_nodes = []
		z_nodes = []
		values_list = []
		#~ print coords_vals # DEBUG
		for item in coords_vals:
			x_nodes.append(item[0])
			y_nodes.append(item[1])
			z_nodes.append(item[2])
			values_list = np.append(values_list, coords_vals[item]) # What is the difference between the old and new values list
		#~ print values_list
		#~ xs = x_nodes
		#~ ys = y_nodes
		#~ zs = z_nodes
		#~ values_list = np.array(values_list)

		#~ lat = self.lat_len
		#~ lon = self.lon_len
		
		#~ ti_time = np.arange(0, self.time_len)
		#~ ti_lat = np.arange(0, self.lat_len)
		#~ ti_lon = np.arange(0, self.lon_len)

		#~ lati, loni = np.meshgrid(ti_lat, ti_lon)

		#~ s = len(x_nodes)
		
		### TRY THIS +++++++++++++++++++++++++++++++++++
		### Trying out 3d RBF (insert into ga_class)
		### meshgrid 3d
		xxx, yyy, zzz = np.lib.index_tricks.mgrid[0:10, 0:40, 0:180]


		# use RBF
		rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)
		#~ rbf = Rbf(zs, ys, values_list, function='gaussian', epsilon=4)
		ZI = rbf(xxx, yyy, zzz)
		#~ ZI = rbf(loni, lati)
		#~ ZI = np.swapaxes(ZI, 1, 0)
		ZI = ZI**mask[0, 0:10, :, :]
		### +++++++++++++++++++++++++++++++++++++++++
		
		
		# use RBF
		#~ rbf = Rbf(x_nodes, y_nodes, z_nodes, values_list, function='gaussian', epsilon=4)
		#~ ZI = rbf(ti_time, ti_lat, ti_lon) # seems to work!
		#~ ZI = ZI*np.swapaxes(mask[0, 0, :, :], 1, 0)
		#~ self.array = ma.masked_values(self.array, 1e20)
		
		
		self.array = self.array*mask[0, :, :, :]
		GI = self.array[:, :, :]
		#~ GI = np.reshape(GI, (1, 40, 180))
		#~ GI = np.swapaxes(GI, 1, 0)
		self.diff_array = (GI - ZI)**2
		fitness = np.sqrt(np.sum((GI - ZI)**2))
		#~ print np.mean(GI), np.mean(ZI), np.std(GI), np.std(ZI)
		#~ print "The fitness of the solution is: %g" % fitness
		self.population[index]['fitness'] = fitness
		
	def add_fitness_key(self):
		'''
		Goes through each chromosome in the population and adds 
		a fitness key and a value for that key.
		Used in sort_pop method.
		
		'''
		for i in range(len(self.population)):
			coords_vals = self.calc_chrom_rbf_part_1(i)
			self.calc_chrom_rbf_part_2(i, coords_vals)
			# print i # debug

		
		
	def sort_pop(self):
		'''
		Uses the fitness key of the chromosomes
		in the poputlation to sort the population form 
		fit to least fit.
		'''
		self.population = sorted(self.population, key=itemgetter('fitness')) 
		self.fittest_fitness = self.population[0]['fitness']
		self.least_fitness = self.population[-1]['fitness']
		print "The fittest individual has a fitness of %g. The least fit individual has a fitness of %g" % (self.fittest_fitness, self.least_fitness)
		print '*'*79
		
	#~ def roulette_sum(self):
		#~ '''
		#~ Use for Roulette Wheel Selection.
		#~ '''
		
		#~ sum_fitness = 0
		#~ for i in range(len(self.population)):
			#~ sum_fitness = sum_fitness + self.population[i]['fitness']
		#~ self.sum_fitness = sum_fitness

	#~ def roulette_rand(self):
		#~ return random.random()*self.sum_fitness

	#~ def roulette_sel(self):
		#~ self.roulette_sum()
		#~ random_number = self.roulette_rand()
		#~ accrued_sum = 0
		#~ for i in range(len(self.population)):
			#~ accrued_sum = accrued_sum + self.population[i]['fitness']
			#~ if accrued_sum >= random_number:
				#~ parent = self.population[i]
				#~ break

		#~ return self.population[i]['chrom_list']

	def tourn_sel(self):
		'''
		Chooses two random individuals from the
		population.  
		The individual with the lowest value for the fitness key
		becomes a parent for the next generation.
		'''
		x = random.randint(0, 19)
		player1 = self.population[x]
		y = random.randint(0, 19)
		player2 = self.population[y]
		if player1['fitness'] <= player2['fitness']:
			parent = player1['chrom_list']
		else:
			parent = player2['chrom_list']
		return parent
		print x, y

	#~ def select_parent_from_roulette(self):
		#~ return self.roulette_sel()
		
	def select_parent_from_tournament(self):
		'''
		Selects a parent using tournament selection.
		'''
		return self.tourn_sel()
		
	def crossover(self):
		'''
		Selects two parents from the population
		and then mutates the parents
		and then creates two children from 
		the two parents using crossover.
		'''
		crossover_point = random.randint(0, len(self.population[0]['chrom_list']))
		# print 'crossover on %i' % crossover_point
		parent1 = self.select_parent_from_tournament()
		parent2 = self.select_parent_from_tournament()
		parent1 = self.mutate(parent1)
		parent2 = self.mutate(parent2)
		child1 = parent1[:crossover_point] + parent2[crossover_point:]
		child2 = parent1[crossover_point:] + parent2[:crossover_point]
		return child1, child2
		
	def mutate(self, chromosome):
		# self.chromosome = []
		for i in range(len(chromosome)):
			if random.random() < self.mutation_rate:
				print 'mutation on %i' % i
				if chromosome[i] =='0':
					chromosome[i] = '1'
				else:
					chromosome[i] = '0'
		return chromosome
		
	def make_new_pop(self):
		'''
		Uses the crossover function ten times
		in order to get twenty children
		and make a new population.
		'''
		self.new_pop = []
		for i in range(10):
			dictionary1= {}
			dictionary2 = {}
			dictionary1['chrom_list'], dictionary2['chrom_list'] = \
			self.crossover()
			self.new_pop = np.append(self.new_pop, [dictionary1, dictionary2])
			
	def elitism(self):
		'''
		Selects a random individual from the new population
		and replaces it with the fittest.
		'''
		r = random.randint(0, 19)
		self.new_pop[r] = self.population[0]
		### Double elitism
		r = random.randint(0, 19)
		self.new_pop[r] = self.population[1]
		
### --------- OLD CODE ---------- ###
'''
	def calc_chrom_1(self, index):
		chromosome = self.population[index]['chrom_list']
		gene_stepper = 0
		values_list = []
		coord_list = []
		gene_length = self.string_length # hard code
		# print gene_length
		while gene_stepper < len(chromosome) - gene_length:
			gene_list = chromosome[gene_stepper:gene_stepper+20]
			current_gene = ''
			bit_stepper = 0
			while bit_stepper < gene_length:
				current_gene = current_gene + gene_list[bit_stepper]
				bit_stepper += 1
			values_list = np.append(values_list, self.gene_map[current_gene]['value'])
			gene_stepper += gene_length
		fitness_mean = np.abs(ma.mean(ma.masked_values(values_list, 1e+20)) - ma.mean(self.array))
		fitness_stdev = np.abs(ma.std(ma.masked_values(values_list, 1e20)) - ma.std(self.array))
		fitness_masked_count = ma.count_masked(ma.masked_values(values_list, 1e20))
		fitness_1 = 100-((fitness_mean+fitness_stdev)*fitness_masked_count/10)
		fitness_2 = fitness_mean+fitness_stdev #+np.float(fitness_masked_count/1000)
		# print "The fitness of the chromosome/solution is %g." % np.add(fitness_mean, fitness_stdev)
		#if (fitness_mean<self.fitness_mean):
		#~ self.fitness_mean = fitness_mean
		#if (fitness_stdev < self.fitness_stdev):
		#~ self.fitness_stdev = fitness_stdev
		#if (fitness_masked_count <self.fitness_masked_count ):
		#~ self.fitness_masked_count = fitness_masked_count
		#~ self.population[index]['fitness'] = fitness_2
		#~ self.population[index]['fitness_mean'] = self.fitness_mean
		#~ self.population[index]['fitness_stdev'] = self.fitness_stdev
		#~ self.population[index]['fitness_masked_count'] = self.fitness_masked_count
		# print len(values_list) # debug
		
	def calc_chrom_2(self, index):
		dictionary = self.gene_map
		chromosome = self.population[index]['chrom_list']
		gene_stepper = 0
		values_list = []
		coord_list = [(0, 0, 1)]
		gene_length = self.string_length-1 # hard code
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
		#~ print coord_list, values_list
		###
		fitness_mean = np.abs(ma.mean(ma.masked_values(values_list, 1e06)) - ma.mean(self.array))
		fitness_stdev = np.abs(ma.std(ma.masked_values(values_list, 1e06)) - ma.std(self.array))
		fitness_masked_count = ma.count_masked(ma.masked_values(values_list, 1e06))
		fitness_1 = 100-((fitness_mean+fitness_stdev)*fitness_masked_count/10)
		fitness_2 = fitness_mean+fitness_stdev #+np.float(fitness_masked_count/1000)
		# print "The fitness of the chromosome/solution is %g." % np.add(fitness_mean, fitness_stdev)
		#if (fitness_mean<self.fitness_mean):
		self.fitness_mean = fitness_mean
		#if (fitness_stdev < self.fitness_stdev):
		self.fitness_stdev = fitness_stdev
		#if (fitness_masked_count <self.fitness_masked_count ):
		self.fitness_masked_count = fitness_masked_count
		#self.population[index]['fitness'] = fitness_2
		self.population[index]['fitness_mean'] = self.fitness_mean
		self.population[index]['fitness_stdev'] = self.fitness_stdev
		self.population[index]['fitness_masked_count'] = self.fitness_masked_count
		# print len(values_list) # debug
		###
		#coord_list, values_list = calc_chrom(solution, genetic.gene_map)

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
			
			
		coord_list = np.unique(coords_vals)
		values_list = coords_vals.values()

		xs = x_nodes
		ys = y_nodes
		zs = z_nodes
		values_list = np.array(values_list)

		lat = self.lat_len
		lon = self.lon_len

		ti_lat = np.arange(0, lat)
		ti_lon = np.arange(0, lon)

		lati, loni = np.meshgrid(ti_lat, ti_lon)

		s = len(xs)

		# use RBF
		rbf = Rbf(zs, ys, values_list, epsilon=2)
		ZI = rbf(loni, lati)
		self.array = ma.masked_values(self.array, 1e20)
		GI = ma.mean(self.array[:, :, :], axis=1)
		#~ GI = np.reshape(GI, (1, 40, 180))
		GI = np.swapaxes(GI, 1, 0)
		self.diff_array = (GI - ZI)**2
		fitness = 500 - np.sqrt(np.sum((GI - ZI)**2))
		#print "The fitness of the solution is: %g" % fitness
		self.population[index]['fitness'] = fitness
 
	def calc_chrom_rbf_part_1(self, index):

		#~ Calculates the fitness of a chromosome according to 
		#~ the square of the difference of an RBF interpolation of
		#~ the sampled points from the chromosome and the
		#~ model data.
		#~ Used in add_fitness_key method.

		

		
		dictionary = self.gene_map
		chromosome = self.population[index]['chrom_list']
		gene_stepper = 0
		values_list = []
		coord_list = [(0, 0, 1)]
		gene_length = self.string_length # hard code ? AGAIN -1 ???
		while gene_stepper < len(chromosome) - gene_length: # The gene stepper goes from gene to gene o chromosome
			gene_list = chromosome[gene_stepper:gene_stepper+20] # chooses a gene
			current_gene = '' # initialises current gene
			bit_stepper = 0 # counter for the bits
			while bit_stepper < gene_length:
				current_gene = current_gene + gene_list[bit_stepper] #constructs the current gene
				bit_stepper += 1
			values_list = np.append(values_list, dictionary[current_gene]['value']) #appends the value of the current gene to a list
			coord_tuple = dictionary[current_gene]['coordinate'] # assigns  current coordinate to a variable
			coord_list.append(coord_tuple) # appends the coordinate of the current gene to the list
			gene_stepper += gene_length # goes to the next gene in the loop
		coord_list.pop(0) # pops out the original coord
		
		coords_vals = dict(zip(coord_list, values_list))
		return coords_vals
		
	def 	calc_chrom_rbf_part_2(self, index, coords_vals):
		
		from Scientific.IO.NetCDF import NetCDFFile
		
		### get landmask
		nc = NetCDFFile('/home/nicholas/data/ORCA2_landmask.nc','r')
		mask = ma.masked_values(nc.variables['MASK'][:, :, :40, :180], -9.99999979e+33)
		nc.close()
		
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
			values_list.append(coords_vals[item]) # What is the difference between the old and new values list

		xs = x_nodes
		ys = y_nodes
		zs = z_nodes
		values_list = np.array(values_list)

		lat = self.lat_len
		lon = self.lon_len

		ti_lat = np.arange(0, lat)
		ti_lon = np.arange(0, lon)

		lati, loni = np.meshgrid(ti_lat, ti_lon)

		s = len(xs)

		# use RBF
		rbf = Rbf(ys, zs, values_list, function='gaussian', epsilon=4)
		ZI = rbf(lati, loni) # seems to work!
		#~ ZI = ZI*np.swapaxes(mask[0, 0, :, :], 1, 0)
		#~ self.array = ma.masked_values(self.array, 1e20)
		self.array = self.array*mask[0, :, :, :]
		GI = self.array[0, :, :]
		#~ GI = np.reshape(GI, (1, 40, 180))
		GI = np.swapaxes(GI, 1, 0)
		self.diff_array = (GI - ZI)**2
		fitness = np.sqrt(np.sum((GI - ZI)**2))
		#~ print np.mean(GI), np.mean(ZI), np.std(GI), np.std(ZI)
		#~ print "The fitness of the solution is: %g" % fitness
		self.population[index]['fitness'] = fitness
'''
