#! usr/bin/python

# roulette_selection.py

'''This module has a function which will take a sorted list of dictionaries. 
The total fitness will be calculated. A random number will be generated. The 
fitnesses will  be summed from 0-n until the sum is greater than r. The 
chromosome last added to the sum will be selected as a parent.

''' 

import random

def roulette_sum(sorted_list):
    sum_fitness = 0
    for i in range(len(sorted_list)):
        sum_fitness = sum_fitness + sorted_list[i]['fitness']
    return sum_fitness

def roulette_rand(sum_fitness):
    return random.random()*sum_fitness

def roulette_sel(sorted_list):
    sum_fitness = roulette_sum(sorted_list)
    random_number = roulette_rand(sum_fitness)
    accrued_sum = 0
    for i in range(len(sorted_list)):
        accrued_sum = accrued_sum + sorted_list[i]['fitness']
        # print accrued_sum, i #debug
        if accrued_sum >= random_number:
            parent = sorted_list[i]
            # print i, sorted_list[i]['fitness'], accrued_sum
            break
    print "random number is %g and sum fitness is %g" % (random_number, sum_fitness) #debug
    print "individual/parent fitness is %g" % (sorted_list[i]['fitness'])
    return parent


