# make_population.py

'''A tool to make a population of chromosomes
that correspond to the gene_map created by the
make_gene_map module.

'''

import random
import numpy as np

class InitPopulation(object):

	def __init__(self):
		self.pop_list = [] # was {}
		self.pop_size = 20
		self.chromosome_size = 200 # hard code
		self.gene_size = 0

 
	def get_gene_size(self, number):
		self.gene_size = number 
		# print self.gene_size # debug

	def make_chromosome(self):
		self.new_chromosome = []
		for i in range(self.gene_size*self.chromosome_size):
			self.new_chromosome.append(str(random.randint(0,1)))
		# print self.new_chromosome # debug
		return self.new_chromosome

	def set_up_new_pop(self):
		list = [self.make_chromosome() for i in range(self.pop_size)]
		for i in range(self.pop_size):
			dictionary = {}
			dictionary['individual_no'] = i
			dictionary['chrom_list'] = list[i]
			self.pop_list = np.append(self.pop_list, dictionary)
