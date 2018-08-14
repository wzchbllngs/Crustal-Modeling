#"Your mind will answer most questions if you learn to relax and wait for the answer."
#  -- Burroughs
#
# v10.0:: Created August 16, 2016.
#
# This script operates slightly differently from all the other scripts written
#	up until now--this one requires an input file for all data being fed into
#	the script. 
#
# Tags:
#	- hole :: the site to be examined (must be included in the data file's name).
#	- dataset :: the specific dataset to be looked at (must be included in the 
#		data file's name).
#	- datapath :: the path to get the file from (excludes the file name). It m
#	- outputpath :: the path the output file will show up in (excludes the file
#		name).
#	- outputname (optional) :: the name of the output files. If not specified,
#		it is automatically set to whatever the input file was named. The
#		extensions for the two files are ".output" (for the actual results) and
#		".plot" (for the data needed to plot it).
#	- tide :: the tidal constituent. As many constituents as one wants can be 
#		included; each one must get its own "tide = " tag. There must be at least
#		one constituent in the input file.
#	- gammastep :: the step size for the loading efficiency (i.e. gamma) when 
#		doing the brute force calculations. Gamma is bound by 0 and 1.
#	- gammafloor (optional) :: the minimum value for the loading efficiency. If
#		not specified it will initialize to the gammastep value given.
#	- gammacap (optional) :: the maximum value for the loading efficiency. If not
#		specified it will initialize to 1.0.
#	- etastep :: the step size for the hydraulic diffusivity (i.e. eta) when doing
#		brute force calculations. It must be positive/non-zero, and smaller than
#		the etacap (see below).
#	- etafloor (optional) :: the minimum value for the hydraulic diffusivity. If
#		not specified it will initialize to the etastep value given.
#	- etacap (optional) :: the maximum value for the hydraulic diffusivity. If not
#		specified it will initialize to 3000.0 square meters per second.
#	- errorcap (optional) :: the maximum value allowed for the error values calculated
#		during the brute force calculations. The normal range of the error is
#		between 0. and 100. %. If not specified, the errorcap will initialize
#		to 50.
#	- Source <number> :: The parameters for a given source. There must be at least
#		one source for the script to run. Two numbers are given for a source, to
#		be separated by a space. The first number is the distance (in meters) of
#		the source to the CORK. The second number is the loading efficiency at the
#		interface between the sediment and the basement; if there is no sediment
#		the second number is not necessary (it will initialize to 1). Various
#		possible distances and interface gammas can be tested; for each source
#		each possibility must be specified. For example, suppose there are three
#		sources, of which the third is sedimented. We can included a maximum and
#		minimum distance for each source to be tried as such, with various gammas
#		included for the third source:
#		Source 1 = 10000;
#		Source 1 = 12000;
#		Source 2 = 6000;
#		Source 2 = 7000;
#		Source 3 = 6500 .71;
#		Source 3 = 7500 .71;
#		Source 3 = 6500 .63;
#		Source 3 = 7500 .63;
#
#	After each tag, an equal sign must be placed in front, with a space separating
#		everything. Each line must end with a semi-colon. Comments CAN be thrown in;
#		these require a hashtag before the line. For example:
#--------------------------------------
#	# Here is a comment in the input file!
#	hole = 1024C;
#	dataset = 1997;
#	datapath = ~/location/hereitis/;
#	outputpath = ~/location/whereitgoes/;
#	tide = M2;
#	tide = S2;
#	tide = O1;
#	gammastep = 0.05;
#	etastep = 10;
#	errorcap = 20;
#	Source 1 = 10000;
#	Source 1 = 12000;
#	Source 2 = 6000;
#	Source 2 = 7000;
#	Source 3 = 6500 .71;
#	Source 3 = 7500 .71;
#	Source 3 = 6500 .63;
#	Source 3 = 7500 .63;
#--------------------------------------
#
# To separate it from past versions of the script, I am starting the count at
#	version 10.
# For the data files, the following info must be included in the file name:
#	- hole
#	- dataset
#	- borehole or seafloor (whether it is recording basement or seafloor data)
#	- either a .txt, .dat, or .data extension.
#
#==============================================================================
# v10.1:: Created September 14, 2016.
#
# This version is basically the same as 10.0. The main differences are:
#	1. The Errror() function will no longer be a function, but will instead
#		be put directly into the script. This is done to help speed things up.
#	2. The pressure and phase errors will be included as well. They will be
#		included in the output file.
#	3. The phase error is calculated slightly differently. In all programs up to
#		this point, the error was found by comparing a calculated phase (from
#		arctangent) to the original phase. Because the arctangent restricts the
#		calculated value to quadrants 1 and 4, a larger difference than what actually
#		is will be calculated if, for example, the actual phase is in quadrants
#		2 or 3. By taking the tangent then arctangent of the actual phase, the
#		value to be compared to for the error will be restricted to quadrants 1 and
#		4 as well.
#	4. There was actually a mistake in the 10.0 Errror() function. the gamma0s list
#		was not included as part of the calculation.
#
#==============================================================================
# v10.2:: Created October 11, 2016.
#
# This version is similar to 10.1. The differences/improvements are:
#	1. The pressure and phase errors are now calculated after the minimum error loading
#		efficiency and hydraulic diffusivity are found, rather than dealing with
#		sorting. This is done in part because I'm not sure it is sorting those errors
#		properly.
#	2. The iteration tuples/lists will now be generated once, and put into a list that will
#		then be looped over (in a for loop) for each tidal constituent. This list will
#		be sorted.
#		The reason for this is so that each constituent has the same cases. While it
#		isn't particularly important for this script itself, it allows for more options
#		and easier scripting for the plotting scripts 2D_plot-10.#.py, particularly
#		ensuring that each constituent has the same source parameters for each case.
#		It also makes the cases consistent every time, no matter how many times this
#		script is run.
#
#==============================================================================
import numpy as np
import tidedata_extractor as tidex
import what_is_it as wit
from potpourri import does_it_end_in as ends
from potpourri import decimalplace_counter as dpc
from potpourri import minfinder
from random import randrange

import matplotlib.pyplot as plt

import os
import sys

print ' '
#==============================================================================
# Gonna handle the input file here, sort the data, etc.
input_file = sys.argv[1]
file = open(input_file, 'r')
line_list = []
for line in file:
	if (line != '\n') and (line[0] != '#'):
		if ';' not in line:
			print "Check and make sure all valid lines end in a semi-colon (;)."
			exit()
		else:
			line_list.append(line)

source_list = []
tide_list = [] # Check length of the list in order to see if anything there.
iteration_dict = {}
gamma0_dict = {}
z_dict = {} # distance dictionary

datapath_ctr = 0
outputpath_ctr = 0
hole_ctr = 0
dataset_ctr = 0
gammastep_ctr = 0
gammafloor_ctr = 0
gammacap_ctr = 0
etastep_ctr = 0
etafloor_ctr = 0
etacap_ctr = 0
errorcap_ctr = 0
outputname_ctr = 0

# If any of these counters is zero, list after that the thing is missing, and exit.

for line in line_list:
	if "hole = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		hole = splitline2[0]
		hole_ctr += 1
		
	if "dataset = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		dataset = splitline2[0]
		dataset_ctr += 1
	
	if "datapath = ".upper() in line.upper():
		splitline1 = line.split('=')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		datapath = splitline2[0]
		datapath = datapath[1:]
		datapath_ctr += 1
	
	if "outputpath = ".upper() in line.upper():
		splitline1 = line.split('=')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		outputpath = splitline2[0]
		outputpath = outputpath[1:]
		outputpath_ctr += 1
	
	if "outputname = ".upper() in line.upper():
		splitline1 = line.split('=')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		outputname = splitline2[0]
		outputname_ctr += 1
	
	if "tide = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		tide = splitline2[0].upper()
		tide_list.append(tide)
	
	if "gammastep = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		gammastep = float(splitline2[0])
		gammastep_ctr += 1
	
	if "gammafloor = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		gammafloor = float(splitline2[0])
		gammafloor_ctr += 1
	
	if "gammacap = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		gammacap = float(splitline2[0])
		gammacap_ctr += 1
		
	if "etastep = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		etastep = float(splitline2[0])
		etastep_ctr += 1

	if "etafloor = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		etafloor = float(splitline2[0])
		etafloor_ctr += 1
	
	if "etacap = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		etacap = float(splitline2[0])
		etacap_ctr += 1
	
	if "errorcap = ".upper() in line.upper():
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		errorcap = float(splitline2[0])
		errorcap_ctr += 1

	src_num = 1
	while src_num <= 100:
		src_label1 = "Source "+str(src_num)+" ="
		if src_label1.upper() in line.upper():
			splitline1 = line.split(' ')
			if len(splitline1) == 4:
				rawinput = splitline1[-1]
				splitline2 = rawinput.split(";")
				z = float(splitline2[0])
				gam0 = 1.
			elif len(splitline1) == 5:
				z = float(splitline1[3])
				rawinput = splitline1[-1]
				splitline2 = rawinput.split(";")
				gam0 = float(splitline2[0])
				
			if src_num not in source_list:
				source_list.append(src_num)
				gamma0_dict[src_num] = [gam0]
				z_dict[src_num] = [z]
				iteration_dict[src_num] = 1
			else:
				gamma0_dict[src_num].append(gam0)
				z_dict[src_num].append(z)
				iteration_dict[src_num] += 1
				#print iteration_dict
		src_num += 1

if hole_ctr == 0:
	print "Please include a hole (i.e. site) in the input file. (hole = )"
	exit()
if hole_ctr > 1:
	print "Please include only one hole/site. (hole = )"
	exit()
	
if dataset_ctr == 0:
	print "Please include a dataset label. (dataset = )"
	exit()
if dataset_ctr > 1:
	print "Please include only one dataset. (dataset = )"
	exit()

if datapath_ctr == 0:
	print "Please include a path to retrieve data from. (datapath = )"
	exit()
if datapath_ctr > 1:
	print "Please include only one data path. (datapath = )"
	exit()
if outputpath_ctr == 0:
	print "Please include a path to put output file to. (outputpath = )"
	exit()
if outputpath_ctr > 1:
	print "Please include only one output path. (outputpath = )"
	exit()

if outputname_ctr == 0:
	namesplit1 = input_file.split('/') # for mac
	name_of_interest = namesplit1[-1]
	namesplit = name_of_interest.split('.')
	output_file = namesplit[0]+".out"
	plot_file = namesplit[0]+".plot"
	points_file = namesplit[0]+".points"
if outputname_ctr > 1:
	print "Please put only one output file name. (outputname = )"
	exit()

if gammastep_ctr == 0:
	print "Include the step size for the loading efficiency. (gammastep = )"
	exit()
if gammastep_ctr > 1:
	print "Include only one specific step size for the loading efficiency. (gammastep = )"
	exit()
if gammafloor_ctr == 0:
	gammafloor = gammastep
elif gammafloor_ctr > 1:
	print "Include only one minimum value for loading efficiency. (gammafloor = )"
	exit()
if gammacap_ctr == 0:
	gammacap = 1.
elif gammacap_ctr > 1:
	print "Include only one maximum value for loading efficiency. (gammacap = )"
	exit()

if etastep_ctr == 0:
	print "Include the step size for the hydraulic diffusivity. (etastep = )"
	exit()
if etastep_ctr > 1:
	print "Include only one specific step size for the hydraulic diffusivity. (etastep = )"
	exit()
if etafloor_ctr == 0:
	etafloor = etastep
elif etafloor_ctr > 1:
	print "Include only one minimum value for the hydraulic diffusivity. (etafloor = )"
	exit()
if etacap_ctr == 0:
	etacap = 3000.
elif etacap_ctr > 1:
	print "Include only one maximum value for the hydraulic diffusivity. (etacap = )"

if errorcap_ctr == 0:
	errorcap = 50. #percent
if errorcap_ctr > 1:
	print "Please include only one maximum error value, in percent. (errorcap = )"
	exit()

if len(tide_list) == 0:
	print "Please include at least one constituent."
	exit()
if len(source_list) == 0:
	print "Please include at least one source."
	exit()

#print input_file
#print output_file
#print dataset
#print datapath
#print tide_list
#print source_list
#print gamma0_dict
#print z_dict

#==============================================================================
# FUNCTIONS:


#==============================================================================
# Collecting the data files:

directory = os.listdir(datapath)
#print directory

ctr = 0
for thing in directory:
	if ends(".txt",thing) or ends(".dat",thing) or ends(".data",thing):
		if (hole in thing) and (dataset in thing):
			if 'borehole' in thing:
				bfile = thing
				ctr += 1
			if 'seafloor' in thing:
				sfile = thing
				ctr += 1
if ctr != 2:
	print "Check data files in "+datapath
	exit()
#----------------------------------------------------------
# List creation:

gamlist_init = []
etalist_init = []

gam_ctr = gammafloor
while gam_ctr <= gammacap:
	gamlist_init.append(gam_ctr)
	gam_ctr += gammastep
eta_ctr = etafloor
while eta_ctr <= etacap:
	etalist_init.append(eta_ctr)
	eta_ctr += etastep

#print gamlist_init
#print etalist_init


#----------------------------------------------------------
# Collecting the data from the data files. Since the brute force analysis has to
#	be run for each constituent individually, it will have to loop over the list
#	of tidal constituents.

product = 1
for source in source_list:
	product *= iteration_dict[source]
#print product

final_tide_list = []
gam_dict = {}
eta_dict = {}
error_dict = {}
perr_dict = {}
phierr_dict = {}

min_gam_dict = {}
min_eta_dict = {}
min_error_dict = {}

file_out = open(outputpath+output_file,'w') # main output file
file_plot = open(outputpath+plot_file,'w') # these will be used for any plots wanted
file_points = open(outputpath+points_file,'w')	# This is used to plot ALL points, not just a line. 
												#	This helps when a graph is unclear with the .plot file.

file_out.write("# Source info:: 1. distance in meters (m), 2. loading efficiency at source. \n")
file_out.write("# 'Point' info:: ( loading eff, hydraulic diff (m^2/s), pressure err, phase err, err (RMS percent) ) \n")
file_out.write('\n')
file_plot.write("# Source info:: 1. distance in meters (m), 2. loading efficiency at source. \n")
file_plot.write("# 'Point' info:: ( loading eff, hydraulic diff (m^2/s), pressure err, phase err, err (RMS percent) ) \n")
file_plot.write('\n')
file_points.write("# Loading efficiency, hydraulic diffusivity, error \n")
file_points.write('\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This part here creates the combinations of distances to be plugged into
#	the Errror() function. Because I couldn't think of a loop structure to
#	create all the combinations given the knowns, I have it set up so that
#	it will eventually "guess" all the combinations with pseudo-random
#	numbers.

product = 1
for source in source_list:
	product *= iteration_dict[source]
#print product

iteration_list = []

while len(iteration_list) < product:	
	iteration_combo = []
	for source in source_list:
		i = randrange(0,iteration_dict[source])
		iteration_combo.append(i)
			
	if iteration_combo not in iteration_list:
		iteration_list.append(iteration_combo)

iteration_list.sort()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for tide in tide_list:
	bdata = tidex.basic_extract(datapath+bfile, tide)
	sdata = tidex.basic_extract(datapath+sfile, tide)
	#print bdata
	#print sdata
	T = bdata[0]
	p = bdata[1]
	tide_final = bdata[3]
	print tide_final
	final_tide_list.append(tide_final)
	sigma = sdata[1]
	phi = bdata[2] - sdata[2]

	case = 0	# I'm starting on 0 in case it needs to be used as an index. In
				#	the actual iteration 1 will be added so that case = 1 on start.
	for combo in iteration_list:
		percent = int(float(case)/len(iteration_list)*100)
		print "Roughly "+str(percent)+"% done."

		case += 1 #number of combo that we're on
			
		file_out.write("Case "+str(case)+': '+tide_final+' \n')
		file_plot.write("Case "+str(case)+': '+tide_final+' \n')
		file_points.write("Case "+str(case)+': '+tide_final+' \n')
		# The plotting script and any other script will use the case labels
		#	to encapsulate the cases.
		
		# Including the distances and source gammas for each source:
			
		distances = []
		gamma0s = []
		index = 0
		while index < len(source_list):
			iteration_index = int(combo[index])
			z = z_dict[source_list[index]][iteration_index]
			gam0 = gamma0_dict[source_list[index]][iteration_index]
			distances.append(z)
			gamma0s.append(gam0)
			file_out.write("Source "+str(index+1)+":: ")
			file_plot.write("Source "+str(index+1)+":: ")
			file_points.write("Source "+str(index+1)+":: ")
			file_out.write(str(z)+'m '+str(gam0)+' \n')
			file_plot.write(str(z)+'m '+str(gam0)+' \n')
			file_points.write(str(z)+'m '+str(gam0)+' \n')
			index += 1
#		print distances, gamma0s
#		print iteration_trash
		tmp_gamlist = []
		tmp_etalist = []
		tmp_errlist = []
		for gam in gamlist_init:
			for eta in etalist_init:
# Error calculations:

				cos_sum = 0.
				sin_sum = 0.
	
				ijk = 0
				while ijk < len(distances):
					z = distances[ijk]
					LL = np.pi*z/np.sqrt(np.pi*eta*T)
					gam0 = gamma0s[ijk]
					costerm = (gam0-gam)*np.exp(-LL)*np.cos(LL)
					sinterm = (gam0-gam)*np.exp(-LL)*np.sin(LL)

					cos_sum += costerm
					sin_sum += sinterm
					ijk += 1

				cos_calc = gam + cos_sum
				sin_calc = sin_sum
				cos_insitu = p*np.cos(phi)/sigma
				sin_insitu = p*np.sin(phi)/sigma
	
				cos_error = np.abs(cos_calc - cos_insitu)*100./cos_insitu
				sin_error = np.abs(sin_calc - sin_insitu)*100./sin_insitu

# Calculates the "pseudovector" total error between the real and imaginary
#	errors (Basically does a root mean square). This can also be done
#	with the pressure and phase errors, but real/imaginary was done for
#	simplicity.
				errror = np.sqrt((cos_error**2 + sin_error**2)/2.)
				# These are basically two different errors that can be used--either with cos and sin,
				#	or with the actual phase and pressure errors. While both are being calculated,
				#	for now I will only be doing everything with errror only.
					
				if errror < errorcap:
					tmp_gamlist.append(gam)
					tmp_etalist.append(eta)
					tmp_errlist.append(errror)
					file_points.write(str(gam)+','+str(eta)+','+str(errror)+' \n')
		mins_plot = minfinder(tmp_gamlist,tmp_etalist,tmp_errlist)
		mins_out = minfinder(mins_plot[0],mins_plot[1],mins_plot[2])
		#print mins_out
		
		ii = 0
		while ii < len(mins_plot[0]):
			pt_gam = mins_plot[0][ii]
			pt_eta = mins_plot[1][ii]
			pt_err = mins_plot[2][ii]

			cos_sum = 0.
			sin_sum = 0.

			ijk = 0
			while ijk < len(distances):
				z = distances[ijk]
				LL = np.pi*z/np.sqrt(np.pi*pt_eta*T)
				gam0 = gamma0s[ijk]
				costerm = (gam0-pt_gam)*np.exp(-LL)*np.cos(LL)
				sinterm = (gam0-pt_gam)*np.exp(-LL)*np.sin(LL)
				cos_sum += costerm
				sin_sum += sinterm
				ijk += 1
			cos_calc = pt_gam + cos_sum
			sin_calc = sin_sum
		
			ratio_calc = np.sqrt(cos_calc**2 + sin_calc**2)
			ratio_insitu = p/sigma
			phi_calc = np.arctan(sin_calc/cos_calc)
			phi_compare = np.arctan(np.tan(phi)) # This is included to make sure the
												 #	value calculated for phi can be
												 #	compared more fairly to the original,
												 #	as phi_calc will be restricted to 
												 #	quadrants 1 and 4.
												 # By putting the actual phi through the
												 #	same calculation, it will ensure that
												 #	there's no issue with, say, phi being
												 #	in quadrant 2 and phi_calc being in 4,
												 #	giving a larger error than what actually
												 #	may be.
					
			p_error = np.abs(ratio_calc - ratio_insitu)/ratio_insitu * 100.
			phi_error = np.abs((phi_compare - phi_calc)/phi_compare)*100.

			point = str(pt_gam)+','+str(pt_eta)+','+str(p_error)+','+str(phi_error)+','+str(pt_err)+' \n'
			file_plot.write(point)
			ii += 1

		iii = 0
		while iii < len(mins_out[0]):
			pt_gam = mins_out[0][iii]
			pt_eta = mins_out[1][iii]
			pt_err = mins_out[2][iii]
			
			cos_sum = 0.
			sin_sum = 0.

			ijk = 0
			while ijk < len(distances):
				z = distances[ijk]
				LL = np.pi*z/np.sqrt(np.pi*pt_eta*T)
				gam0 = gamma0s[ijk]
				costerm = (gam0-pt_gam)*np.exp(-LL)*np.cos(LL)
				sinterm = (gam0-pt_gam)*np.exp(-LL)*np.sin(LL)
				cos_sum += costerm
				sin_sum += sinterm
				ijk += 1
			cos_calc = pt_gam + cos_sum
			sin_calc = sin_sum
		

		
			ratio_calc = np.sqrt(cos_calc**2 + sin_calc**2)
			ratio_insitu = p/sigma
			phi_calc = np.arctan(sin_calc/cos_calc)
			phi_compare = np.arctan(np.tan(phi)) # This is included to make sure the
												 #	value calculated for phi can be
												 #	compared more fairly to the original,
												 #	as phi_calc will be restricted to 
												 #	quadrants 1 and 4.
												 # By putting the actual phi through the
												 #	same calculation, it will ensure that
												 #	there's no issue with, say, phi being
												 #	in quadrant 2 and phi_calc being in 4,
												 #	giving a larger error than what actually
												 #	may be.
					
			p_error = np.abs(ratio_calc - ratio_insitu)/ratio_insitu * 100.
			phi_error = np.abs((phi_compare - phi_calc)/phi_compare)*100.
			point = str(pt_gam)+','+str(pt_eta)+','+str(p_error)+','+str(phi_error)+','+str(pt_err)+' \n'
			file_out.write(point)
			iii += 1
			
		file_out.write("End " +str(case)+'\n')
		file_out.write('\n')
		file_plot.write("End "+str(case)+'\n')
		file_plot.write('\n')
		file_points.write("End "+str(case)+'\n')
		file_points.write('\n')

		if case == len(iteration_list):
			print "100% done."

file_out.close()
file_plot.close()
print ' '
print "~~~~~~~~~~~~~~~~~~~~"
print ' '




