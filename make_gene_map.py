#make_gene_map.py

'''Make a gene map that maps a binary string
to a coordinate of an array as well as a value 
for that coordinate.'''

import numpy as np
import itertools
import load_cflux_masked
from numpy import ma

class GeneMap(object):

	def __init__(self):
		self.gene_map = {}
		self.array = []
		self.array_size = 0
		self.time_len = 0
		self.lat_len = 0
		self.lon_len = 0
		self.last_valid_binary_string = ''
		self.string_length = 0
		self.count = 0

	def get_array_attributes(self):
		file_name = '/home/nicholas/data/CFLX_2000_2009.nc'
		time_start = 0
		time_end = 730
		lat_start = 0
		lat_end = 40
		lon_start = 0
		lon_end = 182
		masked_value = 1e20
		self.array = load_cflux_masked.load_file(time_end=1, lat_start = 0, lat_end=40, lon_end=180) #!
		self.array_shape = np.shape(self.array)
		self.array_size = np.size(self.array)
		self.time_len = self.array_shape[0]
		self.lat_len = self.array_shape[1]
		self.lon_len = self.array_shape[2]
		self.string_length = len(bin(self.array_size)[2:])
		print self.string_length
		
	def make_gene_map_1(self):
		count = 0
		self.iterator = itertools.product(range(self.time_len), range(self.lat_len),
									range(self.lon_len))
		for x_valid in self.iterator:
			binary_string = bin(count)[2:]
			while len(binary_string) < self.string_length:
			    binary_string = '0' + binary_string
			# self.gene_map[binary_string] = {}
			# self.gene_map[binary_string]['coordinate'] = tuple(x_valid)
			if ma.getdata(self.array[x_valid]) < 1000:
				self.gene_map[binary_string] = {}
				self.gene_map[binary_string]['coordinate'] = tuple(x_valid)
				self.gene_map[binary_string]['value'] = self.array[x_valid]
				#count += 1
				print count, binary_string, self.gene_map[binary_string]['value'] # debug
			else:
				# self.gene_map[binary_string]['value'] = 1E06
				pass
			count += 1
			# print count, binary_string, self.gene_map[binary_string]['value'] # debug
		self.last_valid_binary_string = binary_string
		not_valid_first = eval('0b' + binary_string) + 1
		not_valid_last = eval('0b' + '1'*self.string_length)
		print not_valid_last
		print x_valid
		
		for x_not_valid in range(not_valid_first, not_valid_last+1):
			binary_string = bin(x_not_valid)[2:]
			self.gene_map[binary_string] = {}
			self.gene_map[binary_string]['coordinate'] = (999, 999, 999)
			self.gene_map[binary_string]['value'] = 1e06
			print x_not_valid, binary_string, self.gene_map[binary_string]['value']
			# print x_not_valid, binary_string # debug

	def make_gene_map_2(self):
		self.array = ma.getdata(self.array)
		count = 0
		self.iterator = itertools.product(range(self.time_len), range(self.lat_len),
									range(self.lon_len))
		for x_valid in self.iterator:
			# print x # debug
			binary_string = bin(count)[2:]
			while len(binary_string) < self.string_length-1:
				# print len(binary_string) #debug
				binary_string = '0' + binary_string
			# print binary_string # debug
			self.gene_map[binary_string] = {}
			# print self.array[x_valid]
			if self.array[x_valid] > 10:
				print "masked"
				pass
			else:
				self.gene_map[binary_string]['coordinate'] = tuple(x_valid)
				self.gene_map[binary_string]['value'] = self.array[x_valid]
				count += 1
				print binary_string, tuple(x_valid), self.array[x_valid]
			# print self.count # += 1
		self.last_valid_binary_string = binary_string
		not_valid_first = eval('0b' + binary_string) + 1
		not_valid_last = eval('0b' + '1'*(self.string_length-1)) # added minus one just for nonmasked version
		self.count = count
		print count
		print binary_string
		print len(binary_string)
		print not_valid_first
		print not_valid_last
		print x_valid
		
		for x_not_valid in range(not_valid_first, not_valid_last+1):
			binary_string = bin(x_not_valid)[2:]
			self.gene_map[binary_string] = {}
			self.gene_map[binary_string]['coordinate'] = (999, 999, 999)
			self.gene_map[binary_string]['value'] = 1e06
			#print x_not_valid, binary_string, self.gene_map[binary_string]['value']
			# print x_not_valid, binary_string # debug