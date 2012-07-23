#! usr/bin/python

# ranking.py

'''Changes a dictionary into a list of dictionaries and then sorts the list.

'''

import numpy as np
from operator import itemgetter

def make_list(dictionary):
    x=0
    list_of_dict = []
    for item in dictionary:
        list_of_dict.append(dictionary[item])
        x+=1
    return list_of_dict

def sort_list(a_list):
    newlist = sorted(a_list, key=itemgetter('fitness'), reverse=False)
    return newlist
