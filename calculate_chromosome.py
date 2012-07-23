# calculate_chromosome.py

''' This Module calculates a chromosome from a dictionary.

'''

#import numpy as np

def calc_chrom(chromosome, dictionary):
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
		#values_list = np.append(values_list, dictionary[current_gene]['value'])
		coord_tuple = dictionary[current_gene]['coordinate']
		print coord_tuple
		coord_list.append(coord_tuple)
		gene_stepper += gene_length
	coord_list.pop(0)
	return values_list, coord_list
