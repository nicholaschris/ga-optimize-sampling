# selection_tools.py

'''This module selects chromosomes from the population dictionary 
at random and through crossover and mutation creates children.

'''

import random

class NewPop(object):
	
	def __init__(self):
		self.mutation_rate = 0
		self.old_pop = {}
		self.new_pop = {}
		self.chromosome1 = []
		self.chromosome2 = []
		
	def mutate(self):
		self.chromosome = []
		for i in range(len(chromosome)):
			if random.random() < mutation_rate:
				print 'mutation on %i' % i
				if chromosome[i] =='0':
					chromosome[i] = '1'
				else:
					chromosome[i] = '0'
		return chromosome
		
	def crossover(self,):
		crossover_point = random.randint(0, len(self.chromosome1))
		print 'crossover on %i' % crossover_point
		child1 = self.chromosome1[:crossover_point] + self.chromosome2[crossover_point:]
		child2 = self.chromosome1[crossover_point:] + self.chromosome2[:crossover_point]
		return child1, child2
		
	def select_parents(self):
		random_number_1 = str(random.randint(0, 19))
		random_number_2 = str(random.randint(0, 19))
		self.chromosome1 = self.old_pop['parent_' + random_number_1]['chrom_list']
		self.chromosome2 = self.old_pop['parent_' + random_number_2]['chrom_list']
		return self.chromosome1, self.chromosome2

	def make_new_pop(self):
		self.select_parent()
		for n in range(1,20, 2):
			self.new_pop['child_' + str(n)] = {}
			self.new_pop['child_' + str(n+1)] = {}
			self.new_pop['child_' + str(n)]['chrom_list'], self.new_pop['child_' + str(n+1)]['chrom_list'] = \
			self.crossover()
		return self.new_pop

    def select_parent_pairs(self, pair):
        ind_1 = pair
        ind_2 = pair + 1
        self.chromosome1 = self.old_pop['parent_' + random_number_1]['chrom_list']
		self.chromosome2 = self.old_pop['parent_' + random_number_2]['chrom_list']
		return self.chromosome1, self.chromosome2

	def crossover_pairs(self):
        for i in range(1, 21, 2):
		    self.select_parent_pairs(i)
		    self.new_pop['child_' + str(i)] = {}
		    self.new_pop['child_' + str(i+1)] = {}
		    self.new_pop['child_' + str(i)]['chrom_list'], self.new_pop['child_' + str(i+1)]['chrom_list'] = \
			self.crossover()
		return self.new_pop
		
# tournament_selection.py

    def pick_parent():
	    x= random.randint(0, 19)
	    player1 = population(x)
	    player2 = population(x)
	    if player1['fitness'] >= player2['fitness']:
		    parent = player1
	    else:
		    parent = parent2
	    return parent
