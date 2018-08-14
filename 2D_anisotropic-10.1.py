# Created September 16, 2016.
#
# This script requires an input file to run, the same way as 2D_isotropic-10.
#	Because it solves for the hydraulic diffusivities parallel and perpendicular
#	to the ridge axis, the loading efficiency at the site has to be given.
#	Due to this, a slightly different set of tags is needed/used.
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
#	- tide :: the tidal constituent information. This tag will have two things 
#		included: the actual tidal constituent, and the loading efficiency for
#		that specific tidal constituent at the site/hole. Both must be entered.
#		This will look something like this:
#			tide = M2 .47;
#		The loading efficiency must be between 0 and 1.
#		As many constituents as one wants can be included; each one must get its 
#		own "tide = " tag. There must be at least one constituent in the input file.
#	Because there are two directions being solved for the 2D hydraulic diffusivity,
#		there can be an etastep, etafloor, and etacap tag for EACH direction. X is the
#		direction perpendicular to the ridge axis, Y is the direction parallel to the
#		ridge axis. If a specific set of tags for X or Y (or both) are listed, they will
#		take priority for that direction over the general set of "eta" tags.
#	- etastep :: the step size for the hydraulic diffusivity (i.e. eta) when doing
#		brute force calculations. It must be positive/non-zero, and smaller than
#		the etacap (see below).
#	- etafloor (optional) :: the minimum value for the hydraulic diffusivity . If
#		not specified it will initialize to the etastep value given.
#	- etacap (optional) :: the maximum value for the hydraulic diffusivity. If not
#		specified it will initialize to 3000.0 square meters per second.
#	- etastep_X, etafloor_X, etacap_X :: same as the above three, but they will be
#		applied only to the X direction.
#	- etastep_Y, etafloor_Y, etacap_Y :: same as that for the X direction, only applied
#		to the Y direction.
#	- errorcap (optional) :: the maximum value allowed for the error values calculated
#		during the brute force calculations. The normal range of the error is
#		between 0. and 100. %. If not specified, the errorcap will initialize
#		to 50.
#	- Source <number> :: The parameters for a given source. There must be at least
#		one source for the script to run. Three numbers are given for a source, to be 
#		separated by a space.
#		The first number is the distance (in meters) of the source
#		to the CORK. It must be greater than 0.
#		The second number is the angle between the signal propagation
#		direction and the normal to the ridge axis IN DEGREES (it will be converted
#		to radians in the script). This should be between 0 and 90 degrees.
#		The third number is the loading efficiency at the interface between the 
#		sediment and the basement; if there is no sediment the third number is not necessary
#		(it will initialize to 1).
#		Various possible distances and interface gammas can be tested; for each source
#		each possibility must be specified. For example, suppose there are three
#		sources, of which the third is sedimented. We can included a maximum and
#		minimum distance for each source to be tried as such, with various gammas
#		included for the third source:
#		Source 1 = 10000 90;
#		Source 1 = 12000 90;
#		Source 2 = 6000 0;
#		Source 2 = 7000 10;
#		Source 3 = 6500 0 .71;
#		Source 3 = 7500 90 .71;
#		Source 3 = 6500 75 .63;
#		Source 3 = 7500 75 .63;
#
# Certain extra tags are required for the 2D_isotropic-10 scripts that are not needed/used
#	for this script. These include gammastep, gammafloor, and gammacap. The loading
#	efficiency that 2D_isotropic-10 solves for is a given in 2D_anisotropic-10 in the
#	"tide" tag.
# After each tag, an equal sign must be placed in front, with a space separating
#	everything. Each line must end with a semi-colon. Comments CAN be thrown in; these
#	require a hashtag before the line.
#
# Here is an example of an input file
#--------------------------------------
#	# Here is a comment in the input file!
#	hole = 1024C;
#	dataset = 1997;
#	datapath = ~/location/hereitis/;
#	outputpath = ~/location/whereitgoes/;
#	tide = M2 .49;
#	tide = M2 .57;
#	tide = S2 .51;
#	tide = O1 .56;
#	etastep = 10;
#	errorcap = 20;
#	etastep_X = 20;
#	Source 1 = 10000 90;
#	Source 1 = 12000 90;
#	Source 2 = 6000 0;
#	Source 2 = 7000 10;
#	Source 3 = 6500 0 .71;
#	Source 3 = 7500 90 .71;
#	Source 3 = 6500 75 .63;
#	Source 3 = 7500 75 .63;
#--------------------------------------
#
# To separate it from past versions of the script, I am starting the count at
#	version 10.
# For the data files, the following info must be included in the file name:
#	- hole
#	- dataset
#	- borehole or seafloor (whether it is recording basement or seafloor data)
#	- either a .txt, .dat, or .data extension.
#==============================================================================
# v10.1: Created October 13, 2016.
#
# Differences between 10.0 and 10.1:
#	1. The iteration tuples will now be generated once, and put into a list that will
#		then be looped over (in a for loop) for each tidal constituent. This list will
#		be sorted.
#		The reason for this is so that each constituent has the same cases. While it
#		isn't particularly important for this script itself, it allows for more options
#		and easier scripting for the plotting scripts 2D_plot-10.#.py, particularly
#		ensuring that each constituent has the same source parameters for each case.
#		It also makes the cases consistent every time, no matter how many times this
#		script is run.
#	2. A shortcut has been included so that the full filepath doesn't have to be inputted
#		when starting up the program. If you just enter the file, the default filepath
#		will be checked--if a filepath is included, this is skipped.
#	3. A ".points" file is now outputted, besides the .out and .plot file. The .points file
#		basically includes ALL points below the errorcap.
#	 This is done mainly for visibility--sometimes it is hard to see in the .plot file 
#		exactly where the minimums are due to a lack of points. On the flip side,
#		sometimes there are too many points in the .points file to render


#==============================================================================

import numpy as np
from random import randrange
import os
import sys

from tidedata_extractor import basic_extract
import what_is_it as wit
from potpourri import does_it_end_in as ends
from potpourri import decimalplace_counter as dpc
from potpourri import minfinder

#==============================================================================
print "\n ======== \n"

inputfile = sys.argv[1]
#print inputfile
if "/" not in inputfile:
	default_path = "/Users/admin/desktop/tidal analysis/archives/2d_anisotropic-10/"
	inputfile = default_path + inputfile
file = open(inputfile,'r')

hole_ctr = 0
dataset_ctr = 0
datapath_ctr = 0
outputpath_ctr = 0
outputname_ctr = 0

etastep_ctr = 0
etafloor_ctr = 0
etacap_ctr = 0
etastepX_ctr = 0
etafloorX_ctr = 0
etacapX_ctr = 0
etastepY_ctr = 0
etafloorY_ctr = 0
etacapY_ctr = 0
errorcap_ctr = 0

etacap_default = 3000. # m^2/s


raw_tidelist = []
# To store constituents straight from the input file.
# Check the length of this list --> works as a tide_ctr
raw_tidegams = []
# To store loading efficiencies from each constituent straight from the 
#	input file.

tidelist = []
# To store a clean list of the tides (with proper capitalization). Each constituent
#	is stored only once.

# For sources (copied from 2D_isotropic-10.0.py):
source_list = []
iteration_dict = {}
gamma0_dict = {}
z_dict = {} # distance dictionary
angle_dict = {} # in radians


flags = 0
line_ctr = 0
for line in file:
	LINE = line.upper()
	if '#'==line[0]:
		continue
	if line=='\n':
		#print "ENTER"
		continue

	line_ctr += 1
#	print line_ctr, line
#	print line[-2]

	if ';' not in (line[-2]+line[-1]):
		print "Missing ; on line "+str(line_ctr)
		flags += 1
		continue
	if ' ' in (line[-3]+line[-2]+line[-1]):
		print "Please no white space between the entry and the ; on line "+str(line_ctr)
		flags += 1
		continue
	
	if 'hole = '.upper() in LINE:
		linesplit1 = line.split(' ')
#		print linesplit1
#		print linesplit1[-1]
		linesplit2 = linesplit1[-1].split(';')
		hole = linesplit2[0]
		hole_ctr += 1
	
	if 'dataset = '.upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		dataset = linesplit2[0]
		dataset_ctr += 1
	
	if "datapath = ".upper() in LINE:
		linesplit1 = line.split(' = ')
		linesplit2 = linesplit1[-1].split(';')
		datapath = linesplit2[0]
		#print datapath
		datapath_ctr += 1
	
	if 'outputpath = '.upper() in LINE:
		linesplit1 = line.split(' = ')
		linesplit2 = linesplit1[-1].split(';')
		outputpath = linesplit2[0]
		outputpath_ctr += 1
	
	if 'outputname = '.upper() in LINE:
		linesplit1 = line.split(' = ')
		linesplit2 = linesplit1[-1].split(';')
		outputname = linesplit2[0]
		outputname_ctr += 1
	
	if 'tide = '.upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		constituent = linesplit1[-2]
		gam_tide = float(linesplit2[0])
		raw_tidelist.append(constituent)
		raw_tidegams.append(gam_tide)
	
	if "etastep = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etastep = float(linesplit2[0])
		etastep_ctr += 1
	if "etafloor = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etafloor = float(linesplit2[0])
		etafloor_ctr += 1
	if "etacap = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etacap = float(linesplit2[0])
		etacap_ctr += 1
		
	if "etastep_X = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etastep_X = float(linesplit2[0])
		etastepX_ctr += 1
	if "etafloor_X = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etafloor_X = float(linesplit2[0])
		etafloorX_ctr += 1
	if "etacap_X = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etacap_X = float(linesplit2[0])
		etacapX_ctr += 1

	if "etastep_Y = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etastep_Y = float(linesplit2[0])
		etastepY_ctr += 1
	if "etafloor_Y = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etafloor_Y = float(linesplit2[0])
		etafloorY_ctr += 1
	if "etacap_Y = ".upper() in LINE:
		linesplit1 = line.split(' ')
		linesplit2 = linesplit1[-1].split(';')
		etacap_Y = float(linesplit2[0])
		etacapY_ctr += 1

	if "errorcap = ".upper() in LINE:
		splitline1 = line.split(' ')
		rawinput = splitline1[-1]
		splitline2 = rawinput.split(";")
		errorcap = float(splitline2[0])
		errorcap_ctr += 1

# Source sorting from 2D_isotropic-10.0.py:	
	src_num = 1
	while src_num <= 100:
		src_label1 = "Source "+str(src_num)+" ="
		if src_label1.upper() in LINE:
			#print line
			splitline1 = line.split(' ')
			if len(splitline1) == 5:
				raw_angle = splitline1[-1]
				splitline2 = raw_angle.split(';')
				angle_deg = float(splitline2[0])
				angle_rad = angle_deg*np.pi/180.
				z = float(splitline1[-2])
				gam0 = 1.
				#print gam0

			elif len(splitline1) == 6:
				z = float(splitline1[-3])
				angle_deg = float(splitline1[-2])
				angle_rad = angle_deg*np.pi/180.
				rawinput = splitline1[-1]
				splitline2 = rawinput.split(';')
				gam0 = float(splitline2[0])
				#print gam0

			else:
				print "Check input file, sources."
				exit()
				
			if src_num not in source_list:
				source_list.append(src_num)
				gamma0_dict[src_num] = [gam0]
				z_dict[src_num] = [z]
				angle_dict[src_num] = [angle_rad]
				iteration_dict[src_num] = 1
			else:
				gamma0_dict[src_num].append(gam0)
				z_dict[src_num].append(z)
				angle_dict[src_num].append(angle_rad)
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
	namesplit1 = inputfile.split('/') # for mac
	name_of_interest = namesplit1[-1]
	namesplit = name_of_interest.split('.')
	output_file = namesplit[0]+".out"
	plot_file = namesplit[0]+".plot"
	point_file = namesplit[0]+".points"
if outputname_ctr > 1:
	print "Please put only one output file name. (outputname = )"
	exit()

# FIGURING OUT THE HYRDRAULIC DIFFUSIVITY PARAMETERS:
if etastepX_ctr == 0:
	if etastep_ctr == 0:
		print "Please include the step size for the hydraulic diffusivity. (etastep = )"
		exit()
	elif etastep_ctr == 1:
		etastep_X = etastep
	elif etastep_ctr > 1:
		print "Please include only one stepsize for the hydraulic diffusivity. (etastep = )"
		exit()
elif etastepX_ctr == 1:
	pass # exastep_X was set above when the script parses the input file.
elif etastepX_ctr > 1:
	print "Please include only one stepsize for the perpendicular hydraulic " \
	+"diffusivity  (etastep_X = )"
	exit()

if etafloorX_ctr == 0:
	if etafloor_ctr == 0:
		etafloor_X = etastep
	elif etafloor_ctr == 1:
		etafloor_X = etafloor
	elif etafloor_ctr > 1:
		print "Please include only one minimum value for the hydraulic diffusivity. (etafloor = )"
		exit()
elif etafloorX_ctr == 1:
	pass # etafloor_X was set above when the script parsed the input file.
elif etafloorX_ctr > 1:
	print "Please include only one minimum value for the perpendicular hydraulic " \
	+"diffusivity  (etafloor_X = )"
	exit()

if etacapX_ctr == 0:
	if etacap_ctr == 0:
		etacap_X = etacap_default
	elif etacap_ctr == 1:
		etacap_X = etacap
	elif etacap_ctr > 1:
		print "Please include only one maximum value for the hydraulic diffusivity. (etacap = )"
		exit()
elif etacapX_ctr == 1:
	pass
elif etacapX_ctr > 1:
	print "Please include only one maximum value for the perpendicular hydraulic " \
	+"diffusivity  (etacap_X = )"
	exit()

if etastepY_ctr == 0:
	if etastep_ctr == 0:
		print "Please include the step size for the hydraulic diffusivity. (etastep = )"
		exit()
	elif etastep_ctr == 1:
		etastep_Y = etastep
	elif etastep_ctr > 1:
		print "Please include only one stepsize for the hydraulic diffusivity. (etastep = )"
		exit()
elif etastepY_ctr == 1:
	pass # exastep_Y was set above when the script parses the input file.
elif etastepY_ctr > 1:
	print "Please include only one stepsize for the perpendicular hydraulic " \
	+"diffusivity  (etastep_Y = )"
	exit()

if etafloorY_ctr == 0:
	if etafloor_ctr == 0:
		etafloor_Y = etastep
	elif etafloor_ctr == 1:
		etafloor_Y = etafloor
	elif etafloor_ctr > 1:
		print "Please include only one minimum value for the hydraulic diffusivity. (etafloor = )"
		exit()
elif etafloorY_ctr == 1:
	pass # etafloor_Y was set above when the script parsed the input file.
elif etafloorY_ctr > 1:
	print "Please include only one minimum value for the perpendicular hydraulic " \
	+"diffusivity  (etafloor_Y = )"
	exit()

if etacapY_ctr == 0:
	if etacap_ctr == 0:
		etacap_Y = etacap_default
	elif etacap_ctr == 1:
		etacap_Y = etacap
	elif etacap_ctr > 1:
		print "Please include only one maximum value for the hydraulic diffusivity. (etacap = )"
		exit()
elif etacapY_ctr == 1:
	pass
elif etacapY_ctr > 1:
	print "Please include only one maximum value for the perpendicular hydraulic " \
	+"diffusivity  (etacap_Y = )"
	exit()

if errorcap_ctr == 0:
	errorcap = 20.
elif errorcap_ctr > 1:
	print "Please include only one errorcap, or enter nothing for the default value. (errorcap = )"
	exit()
#----------------------------------------------------------------------
# Next up: starting the iterations.

# Opening up the .plot and .out files:
# For v10.1: a .point file is also being included because, in some cases, the plot from
#	the .plot file isn't enough to see what is going on, especially in 3D.
#	The .points file only needs the points and errors, that's it.
file_out = open(outputpath+output_file, 'w')
file_plot = open(outputpath+plot_file, 'w')
file_points = open(outputpath+point_file, 'w')

# Getting data files:
directory = os.listdir(datapath)

dir_check = 0
for obj in directory:
	if (hole in obj) and (dataset in obj):
		if 'borehole' in obj:
			bfile = obj
			dir_check += 1
		if 'seafloor' in obj:
			sfile = obj
			dir_check += 1

if dir_check != 2:
	print "Check input/data files."
	exit()

# List creation:
etaX_raw = []
etaY_raw = []

etaX = etafloor_X
while etaX <= etacap_X:
	etaX_raw.append(etaX)
	etaX += etastep_X
etaY = etafloor_Y
while etaY <= etacap_Y:
	etaY_raw.append(etaY)
	etaY += etastep_Y


tidelist = []
gamlist = []
t = 0

mid_tidelist = []
# To store a list of tides where each is mentioned only once. Because this won't
#	have the proper capitalization, it isn't the final list (hence the "middle" list).
#	This is done so we can do a case for each individual tide--because each constituent
#	can be mentioned more than once in the input file, the way case was counted in
#	2D_isotropic-10.1 won't work here.

mult_dict = {}
case_ctr = {}
# The "mult_dict" is a dictionary containing how many times the tide constituents
#	appear in the input file (i.e. the "multiplicity" of each constituent).
for constituent in raw_tidelist:
	if constituent.upper() not in mid_tidelist:
		mid_tidelist.append(constituent.upper())
		mult_dict[constituent.upper()] = 1
		case_ctr [constituent.upper()] = 1
	if constituent.upper() in mid_tidelist:
		mult_dict[constituent.upper()] += 1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This part here creates the combinations of distances from the sources.
#	Because I couldn't think of a loop structure to create all the combinations 
#	given the knowns, I have it set up so that it will eventually "guess" all 
#	the combinations with pseudo-random numbers.

product = 1
for source in source_list:
	product *= iteration_dict[source]

iteration_list = []
while len(iteration_list) < product:
	iteration_combo = []
	for source in source_list:
		i = randrange(0,iteration_dict[source])
		iteration_combo.append(i)
	if iteration_combo not in iteration_list:
		iteration_list.append(iteration_combo)

iteration_list.sort()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


while t < len(raw_tidelist):
	tide = raw_tidelist[t]
	gam_site = raw_tidegams[t]
	
	bdata = basic_extract(datapath+bfile, tide)
	sdata = basic_extract(datapath+sfile, tide)
	T = bdata[0]
	p = bdata[1]
	sigma = sdata[1]
	phi = bdata[2] - sdata[2]
	tide_final = bdata[3]
	tidelist.append(tide_final)
	gamlist.append(gam_site)
	
	for combo in iteration_list:
		case = case_ctr[tide.upper()]
					
		file_out.write("Case "+str(case)+': '+tide_final+' \n')
		file_plot.write("Case "+str(case)+': '+tide_final+' \n')
		file_points.write("Case "+str(case)+': '+tide_final+' \n')
		print "Case "+str(case)+': '+tide_final
		# The plotting script and any other script will use the case labels
		#	to encapsulate the cases.
			
		# Including the distances and source gammas for each source:

		distances = []
		angles = []
		gamma0s = []
		index = 0
		while index < len(source_list):
			iteration_index = int(combo[index])
			z = z_dict[source_list[index]][iteration_index]
			angle = angle_dict[source_list[index]][iteration_index]
			gam0 = gamma0_dict[source_list[index]][iteration_index]
			distances.append(z)
			angles.append(angle)
			gamma0s.append(gam0)
			file_out.write("Source "+str(index+1)+":: ")
			file_plot.write("Source "+str(index+1)+":: ")
			file_points.write("Source "+str(index+1)+":: ")
			file_out.write(str(z)+'m '+str(angle*180./np.pi)+'deg '+str(gam0)+' \n')
			file_plot.write(str(z)+'m '+str(angle*180./np.pi)+'deg '+str(gam0)+' \n')
			file_points.write(str(z)+'m '+str(angle*180./np.pi)+'deg '+str(gam0)+' \n')
			index += 1
#		print distances, angles, gamma0s
#		print iteration_trash
			
		tmp_etaXlist = []
		tmp_etaYlist = []
		tmp_errlist = []
		for etaX in etaX_raw:
			for etaY in etaY_raw:
				cos_sum = 0
				sin_sum = 0
				ijk = 0
				while ijk < len(distances):
					z = distances[ijk]
					angle = angles[ijk]
					gam0 = gamma0s[ijk]
						
					eta_eff = (np.cos(angle)**2)*etaX + (np.sin(angle)**2)*etaY
					LL = np.pi*z/np.sqrt(np.pi*eta_eff*T)
					costerm = (gam0 - gam_site)*np.exp(-LL)*np.cos(LL)
					sinterm = (gam0 - gam_site)*np.exp(-LL)*np.sin(LL)
						
					cos_sum += costerm
					sin_sum += sinterm
					ijk += 1
				cos_calc = gam_site + cos_sum
				sin_calc = sin_sum
				cos_insitu = p*np.cos(phi)/sigma
				sin_insitu = p*np.sin(phi)/sigma

				cos_error = np.abs(cos_calc - cos_insitu)*100./cos_insitu
				sin_error = np.abs(sin_calc - sin_insitu)*100./sin_insitu

# Calculates the "pseudovector" total error between the real and imaginary
#	errors (Basically does a root mean square). This can also be done
#	with the pressure and phase errors, but real/imaginary was done for
#	simplicity, and less calculation time/simplicity of calculation.

				errror = np.sqrt((cos_error**2 + sin_error**2)/2.)
					# These are basically two different errors that can be used--either with cos and sin,
					#	or with the actual phase and pressure errors. While both are being calculated,
					#	for now I will only be doing everything with errror only.
					
				if errror < errorcap:
					tmp_etaXlist.append(etaX)
					tmp_etaYlist.append(etaY)
					tmp_errlist.append(errror)
					file_points.writelines(str(etaX)+','+str(etaY)+','+str(errror)+' \n')
		minsXY = minfinder(tmp_etaXlist,tmp_etaYlist,tmp_errlist)
		minsXY_tuples = []
		mins_index = 0
		while mins_index < len(minsXY[0]):
			tuple = (minsXY[0][mins_index],minsXY[1][mins_index],minsXY[2][mins_index])
			minsXY_tuples.append(tuple)
			mins_index += 1
					
		tmp_etaXlist = []
		tmp_etaYlist = []
		tmp_errlist = []
		for etaY in etaY_raw:
			for etaX in etaX_raw:
				cos_sum = 0
				sin_sum = 0
				ijk = 0
				while ijk < len(distances):
					z = distances[ijk]
					angle = angles[ijk]
					gam0 = gamma0s[ijk]
						
					eta_eff = (np.cos(angle)**2)*etaX + (np.sin(angle)**2)*etaY
					LL = np.pi*z/np.sqrt(np.pi*eta_eff*T)
					costerm = (gam0 - gam_site)*np.exp(-LL)*np.cos(LL)
					sinterm = (gam0 - gam_site)*np.exp(-LL)*np.sin(LL)
						
					cos_sum += costerm
					sin_sum += sinterm
					ijk += 1
				cos_calc = gam_site + cos_sum
				sin_calc = sin_sum
				cos_insitu = p*np.cos(phi)/sigma
				sin_insitu = p*np.sin(phi)/sigma

				cos_error = np.abs(cos_calc - cos_insitu)*100./cos_insitu
				sin_error = np.abs(sin_calc - sin_insitu)*100./sin_insitu
				errror = np.sqrt((cos_error**2 + sin_error**2)/2.)

				if errror < errorcap:
					tmp_etaXlist.append(etaX)
					tmp_etaYlist.append(etaY)
					tmp_errlist.append(errror)
					file_points.writelines(str(etaX)+','+str(etaY)+','+str(errror)+' \n')
		minsYX = minfinder(tmp_etaYlist, tmp_etaXlist, tmp_errlist)
		minsYX_tuples = []
		mins_index = 0
		while mins_index < len(minsYX[0]):
			tuple = (minsYX[0][mins_index],minsYX[1][mins_index],minsYX[2][mins_index])
			minsYX_tuples.append(tuple)
			mins_index += 1

#	This is to get the points to be plotted. These go into the .plot file.			
		points = []
		for point in minsXY_tuples:
			if point not in points:
				points.append(point)
		for point in minsYX_tuples:
			if point not in points:
				points.append(point)
			
		for point in points:
			cos_sum = 0
			sin_sum = 0
			ijk = 0
			while ijk < len(distances):
				z = distances[ijk]
				angle = angles[ijk]
				gam0 = gamma0s[ijk]
					
				eta_eff = (np.cos(angle)**2)*point[0] + (np.sin(angle)**2)*point[1]
				LL = np.pi*z/np.sqrt(np.pi*eta_eff*T)
				costerm = (gam0 - gam_site)*np.exp(-LL)*np.cos(LL)
				sinterm = (gam0 - gam_site)*np.exp(-LL)*np.sin(LL)
					
				cos_sum += costerm
				sin_sum += sinterm
				ijk += 1
				
			cos_calc = gam_site + cos_sum
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

			plot_point = str(point[0])+','+str(point[1])+','+str(p_error)+','
			plot_point += str(phi_error)+','+str(point[2])+" \n"
			file_plot.write(plot_point)
				

#	This is to get the local minimums from the plot. This goes into the .out file.
		final_minsXY = minfinder(minsXY[0],minsXY[1],minsXY[2])
		final_minsYX = minfinder(minsYX[0],minsYX[1],minsYX[2])

		final_tuplesXY = []
		final_tuplesYX = []
		tuple_index = 0
		while tuple_index < len(final_minsXY[0]):
			tuple = (final_minsXY[0][tuple_index],final_minsXY[1][tuple_index],final_minsXY[2][tuple_index])
			final_tuplesXY.append(tuple)
			tuple_index += 1
		tuple_index = 0
		while tuple_index < len(final_minsYX[0]):
			tuple = (final_minsYX[0][tuple_index],final_minsYX[1][tuple_index],final_minsYX[2][tuple_index])
			final_tuplesYX.append(tuple)
			tuple_index += 1

		final_mins = []
		for value in final_tuplesXY:
			if value not in final_mins:
				final_mins.append(value)
				#print value
		for value in final_tuplesYX:
			if value not in final_mins:
				final_mins.append(value)
				#print value
			
#	Done the same as for the .plot file.
#		print final_mins
		for minimum in final_mins:
			cos_sum = 0
			sin_sum = 0
			ijk = 0
			while ijk < len(distances):
				z = distances[ijk]
				angle = angles[ijk]
				gam0 = gamma0s[ijk]
				
				eta_eff = (np.cos(angle)**2)*point[0] + (np.sin(angle)**2)*point[1]
				LL = np.pi*z/np.sqrt(np.pi*eta_eff*T)
				costerm = (gam0 - gam_site)*np.exp(-LL)*np.cos(LL)
				sinterm = (gam0 - gam_site)*np.exp(-LL)*np.sin(LL)
				
				cos_sum += costerm
				sin_sum += sinterm
				ijk += 1
			
			cos_calc = gam_site + cos_sum
			sin_calc = sin_sum
			ratio_calc = np.sqrt(cos_calc**2 + sin_calc**2)
			ratio_insitu = p/sigma
			phi_calc = np.arctan(sin_calc/cos_calc)
			phi_compare = np.arctan(np.tan(phi))
			p_error = np.abs(ratio_calc - ratio_insitu)/ratio_insitu * 100.
			phi_error = np.abs((phi_compare - phi_calc)/phi_compare)*100.

			plot_point = str(minimum[0])+','+str(minimum[1])+','+str(p_error)+','
			plot_point += str(phi_error)+','+str(minimum[2])+" \n"
			file_out.write(plot_point)
				
		file_out.write("End "+str(case)+'\n')
		file_out.write('\n')
		file_plot.write("End "+str(case)+'\n')
		file_plot.write('\n')
		file_points.write("End "+str(case)+'\n')
		file_points.write('\n')
		case_ctr[tide.upper()] += 1

	t += 1

file_out.close()
file_plot.close()
file_points.close()
#print raw_tidelist
#print raw_tidegams

print ' '
print '~~~~~~~~~~~~~~~~~~~~'
print ' '
	
