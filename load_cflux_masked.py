# load_cflux_masked.py

'''This module loads the cflux netcdf file
 or any netcdf file that needs to be masked.

'''

from numpy import ma
from Scientific.IO.NetCDF import NetCDFFile

# deafult values
file_name = '/home/nicholas/data/CFLX_2000_2009.nc'
time_start = 0
time_end = 730
lat_start = 0
lat_end = 40
lon_start = 0
lon_end = 182
masked_value = 1e+20


def load_file(file_name = file_name, time_start = time_start, 
		      time_end = time_end, lat_start = lat_start, lat_end = lat_end,
	              lon_start = lon_start, lon_end = lon_end, masked_value = masked_value):
		      
	nc = NetCDFFile(file_name, 'r')
	new_array = nc.variables['Cflx'][time_start:time_end, lat_start:lat_end, 
                                     lon_start:lon_end]
	nc.close()
	new_array = ma.masked_values(new_array, masked_value)
	new_array = new_array*1e08
	return new_array    

if __name__ == '__main__':
	load_file()

