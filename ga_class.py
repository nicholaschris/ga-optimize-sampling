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
		
	def add_fitness_key(self):
		for i in range(len(self.population)):
			self.calc_chrom_2(i)
			# print i # debug

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
		self.fitness_mean = fitness_mean
		#if (fitness_stdev < self.fitness_stdev):
		self.fitness_stdev = fitness_stdev
		#if (fitness_masked_count <self.fitness_masked_count ):
		self.fitness_masked_count = fitness_masked_count
		self.population[index]['fitness'] = fitness_2
		self.population[index]['fitness_mean'] = self.fitness_mean
		self.population[index]['fitness_stdev'] = self.fitness_stdev
		self.population[index]['fitness_masked_count'] = self.fitness_masked_count
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
		#return coord_list, values_list
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
			
			
		'''coord_list = np.unique(coords_vals)
		values_list = coords_vals.values()'''

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
		GI = self.array[0, :, :]
		GI = np.swapaxes(GI, 1, 0)
		self.diff_array = (GI - ZI)**2
		fitness = 500 - np.sqrt(np.sum((GI - ZI)**2))
		#print "The fitness of the solution is: %g" % fitness
		self.population[index]['fitness'] = fitness
		
	def sort_pop(self):
		self.population = sorted(self.population, key=itemgetter('fitness'))
		self.fittest_fitness = self.population[0]['fitness']
		self.least_fitness = self.population[-1]['fitness']
		print "The fittest individual has a fitness of %g. The least fit individual has a fitness of %g" % (self.fittest_fitness, self.least_fitness)
		print "fitness_mean: %g *** fitness_stdev: %g *** self.fitness_masked_count: %g" % (self.population[0]['fitness_mean'],self.population[0]['fitness_stdev'],self.population[0]['fitness_masked_count']) 
		#  self.fitness_mean, self.fitness_stdev, np.float(self.fitness_masked_count))
		#print "The least fit individual has a fitness of %g." % self.least_fitness
		# print '*'*79
		
	def roulette_sum(self):
		sum_fitness = 0
		for i in range(len(self.population)):
			sum_fitness = sum_fitness + self.population[i]['fitness']
		self.sum_fitness = sum_fitness

	def roulette_rand(self):
		return random.random()*self.sum_fitness

	def roulette_sel(self):
		self.roulette_sum()
		random_number = self.roulette_rand()
		accrued_sum = 0
		for i in range(len(self.population)):
			accrued_sum = accrued_sum + self.population[i]['fitness']
			if accrued_sum >= random_number:
				parent = self.population[i]
				break
		# print "random number is %g and sum fitness is %g" % (random_number, self.sum_fitness) #debug
		# print "individual/parent fitness is %g" % (self.population[i]['fitness'])
		return self.population[i]['chrom_list']

	def tourn_sel(self):
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

	def select_parent_from_roulette(self):
		return self.roulette_sel()
		
	def select_parent_from_tournament(self):
		return self.tourn_sel()


	def crossover(self):
		crossover_point = random.randint(0, len(self.population[0]['chrom_list']))
		# print 'crossover on %i' % crossover_point
		parent1 = self.select_parent_from_roulette()
		parent2 = self.select_parent_from_roulette()
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
		self.new_pop = []
		for i in range(10):
			dictionary1= {}
			dictionary2 = {}
			dictionary1['chrom_list'], dictionary2['chrom_list'] = \
			self.crossover()
			self.new_pop = np.append(self.new_pop, [dictionary1, dictionary2])
			
	def elitism(self):
		r = random.randint(0, 19)
		self.new_pop[r] = self.population[-1]
		


 


