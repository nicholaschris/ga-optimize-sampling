#! usr/bin/python
# my_ga_aux.py

'''A script that does all the things required for a GA.'''

import make_gene_map
import make_population
import selection_tools
import calculate_chromosome
import calc_fitness
import ranking
import roulette_selection
import random
from numpy import ma

x=0
y=0

### Make the gene-map
my_gene_map = make_gene_map.GeneMap()
my_gene_map.get_array_attributes()
my_gene_map.make_gene_map()

bin_str_len = my_gene_map.string_length
population_array = my_gene_map.array

###Construct the initial population
my_pop = make_population.InitPopulation()
my_pop.get_gene_size(bin_str_len)
my_pop.set_up_new_pop()

old_pop_dict = my_pop.pop_dict

print bin_str_len

sum_x = 0
for item in old_pop_dict:
	print '*'*79
	print item
	values_sample = calculate_chromosome.calc_chrom(old_pop_dict[item]['chrom_list'], my_gene_map.gene_map)[0]
	coords_sample = calculate_chromosome.calc_chrom(old_pop_dict[item]['chrom_list'], my_gene_map.gene_map)[1]
	x = calc_fitness.calc_fitness(values_sample, population_array)
	old_pop_dict[item]['fitness'] = x
	sum_x = sum_x+x
print '-'*79

# rank population
# selection
# select two parents, allow crossover and mutation and get two children
# put them in dictionary
# repeat until dictionary is full

new_pop = selection_tools.NewPop()
new_pop.old_pop = my_pop.pop_dict

new_pop_dict = new_pop.make_new_pop()

sum_y=0
for item in new_pop_dict:
	print '*'*79
	print item
	values_sample = calculate_chromosome.calc_chrom(new_pop_dict[item]['chrom_list'], my_gene_map.gene_map)[0]
	coords_sample = calculate_chromosome.calc_chrom(new_pop_dict[item]['chrom_list'], my_gene_map.gene_map)[1]
	# calc_fitness.calc_fitness(values_sample, population_array)
	y = calc_fitness.calc_fitness(values_sample, population_array)
	new_pop_dict[item]['fitness'] = y
	sum_y=sum_y+y
	
print '-'*79
	
print "The mean value of the array is %g." % ma.mean(my_gene_map.array)
print "The stdev of the array is %g." % ma.std(my_gene_map.array)
print x, y

old_pop_list = ranking.make_list(old_pop_dict)
new_pop_list = ranking.make_list(new_pop_dict)

old_pop_list_sorted = ranking.sort_list(old_pop_list)
new_pop_list_sorted = ranking.sort_list(new_pop_list)

parent1 = roulette_selection.roulette_sel(old_pop_list_sorted)
#print "The fitness of parent1 is: %g" % parent1['fitness']

new_gen = {}
for i in range(20):
    new_gen['individual_' + str(i)] = roulette_selection.roulette_sel(old_pop_list_sorted)
    print new_gen['individual_' + str(i)]['fitness'] # debug

gen_count = 0
while gen_count < 100:
    rank_list = ranking.make_list(new_gen)
    sorted_list = ranking.sort_list(rank_list)
    for i in range(1, 21):
        new_gen['individual_' + str(i)] = roulette_selection.roulette_sel(sorted_list)
    print new_gen['individual_0']['fitness'] # debug
    gen_count += 1
