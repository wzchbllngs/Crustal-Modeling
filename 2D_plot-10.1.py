# Created September 15, 2016.
#
# This script is meant to plot data outputted by the 2D_isotropic-10 and
#	2D_anisotropic-10 python scripts. It takes the .plot and .out files, and
#	plots them, using the error as the y-axis, and the gamma or eta values as
#	the x-axis values.
#
# It only takes three values from those given in the .out and .plot files: it uses
#	the first two values as the x-values (for the first and second plots,
#	respectively). The last value is used as the y-value for both plots.
#	Currently this last value is set to be the error for the real and imaginary
#	parts.
#
# This is a basic plotter. Next versions will allow for more flexibility.
#
# On the prompt line: python 2D_plot-10.0.py <file_path_label> <file_title (exclude extension)>
#
#==============================================================================
# v10.1: Created October 13, 2016.
#
# The update allows you to plot either all the constituents for a given case,
#	or plot all the cases for a given constituent. The minimums for each of these
#	will not be plotted.
# Multiple constituents and cases can be entered as well. Minimum points will be
#	included only for the case of one case.
#
#==============================================================================
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np

#==============================================================================
print '\n ======== \n'

# Just checking to make sure all arguments are included on the command line.
#print len(sys.argv)
if len(sys.argv) == 1:
	print "Please include a file path label and a file name, or a full path to a file."
	exit()
elif len(sys.argv) ==2:
	if len(sys.argv[1])!=1:
		print "Please either include the full path of a file, or a file path label and a file name."
		exit()
elif len(sys.argv)==4:
	print "There are too many arguments in the command line."
	print "Please either include the full path of a file, or a file path label and a file name."
	exit()

# Getting constituents and cases to plot.
tide_raw = []
case_raw = []

print "Enter each constituent of interest, one at a time."
print "Enter 'D' when you are done."
print "Enter 'all' if you would like all constituents included in the input file."
while True:
	tide_input = raw_input(">  ").upper()
	if (tide_input == 'D'):
		break
	if tide_input == 'ALL':
		if len(sys.argv[1])==1:
			break
		else:
			print "To use 'all' option you must change default paths in the program."
			exit()
	tide_raw.append(tide_input)

print ' '
print "Enter each case of interest, one at a time."
print "Enter 'D' when you are done."
print "Enter 'all' if you would like all cases from the .out file."
while True:
	case_input = raw_input(">  ").upper()
	if (case_input == 'D') or (case_input == 'ALL'):
		break
	case_raw.append(case_input)

# Details to get the correct files:
file_path_label = sys.argv[1].upper()
file_title = sys.argv[2]

if file_path_label.upper() == 'A':
	file_path = "/Users/admin/desktop/tidal analysis/archives/2d_anisotropic-10/"
elif file_path_label.upper() == 'I':
	file_path = "/Users/admin/desktop/tidal analysis/archives/2d_isotropic-10/"


input_file = file_path+file_title+'.input'
output_file = file_path+file_title+'.out'
plot_file = file_path+file_title+'.plot'
points_file = file_path+file_title+'.points'

#==============================================================================
# FUNCTIONS:
def get_points(file, tide, case):
	open_file = open(file,'r')
	
	points = []
	
	switch = 0
	for line in open_file:
		LINE = line.upper()
		if ('#' in LINE) or (LINE == '\n'):
			continue
		if ("End "+case).upper() in LINE:
			switch = 0
		if switch == 1:
			if ("Source".upper() in LINE) or ("Loading efficiency".upper() in LINE):
				continue
			else:
				points.append(line)
				#print line
		if (tide.upper() in LINE) and (('Case '+case).upper()+':' in LINE):
			switch = 1

	open_file.close()
	return(points)

def point_to_list(points):
	x1 = []
	x2 = []
	y = []
	for point in points:
		#print point
		point_split = point.split(',')
		x1_value = float(point_split[0])
		x2_value = float(point_split[1])
		err_split = point_split[-1].split(' ')
		y_value = float(err_split[0])
		
		x1.append(x1_value)
		x2.append(x2_value)
		y.append(y_value)
		
	return(x1, x2, y)


#==============================================================================

# Getting final list of tides to work with:
tide_list = []
if tide_input == 'ALL':
#	if case_input == 'ALL':
#		print "Currently cannot have all cases if including all tide constituents."
#		exit()

	if file_path_label.upper() == 'I':
		opened_input = open(input_file, 'r')
		for line in opened_input:
			if line[0] != '#':
				if 'tide' in line:
					linesplit = line.split(' ')
					tidesplit = linesplit[2].split(';')
					tide = tidesplit[0]
					tide_list.append(tide)
		opened_input.close()

	elif file_path_label.upper() == 'A':
		opened_input = open(input_file, 'r')
		for line in opened_input:
			if line[0] != '#':
				if 'tide' in line:
					linesplit = line.split(' ')
					tide = linesplit[2]
					if tide not in tide_list:
						tide_list.append(tide)
		opened_input.close()
	print tide_list

else:
	tide_list = tide_raw

# Getting final list of cases to work with:
case_list = []
if case_input == 'ALL':	
	case_ctr = 0
	opened_out = open(output_file, 'r')
	line_ctr = 0
	
	for line in opened_out:
		line_ctr += 1
		if line[0] != '#':
			if 'Case' in line:
				linesplit = line.split(' ')
				casesplit = linesplit[1].split(':')
				case = int(casesplit[0])
				
				if case > case_ctr:
					case_ctr = case
				elif case == case_ctr:
					print "We've got a problem--a case repeated within the output file!"
					print "Line:" +str(line_ctr)
					print line
				elif case < case_ctr:
					break
	i = 1
	while i <= case_ctr:
		case_list.append(str(i))
		i += 1
	opened_out.close()
	print "Number of cases: "+str(case_ctr)
else:
	case_list = case_raw


#PLOTTING!
if len(case_list) > 1:
	choice1 = raw_input("(M)ins or p(O)ints?  ").upper()
elif len(case_list) == 1:
	choice1 = 'O'

fig,(a0,a1) = plt.subplots(1,2, sharey = True)

if file_path_label == 'I':
	a0.set_xlabel("Loading efficiency")
	a0.set_ylabel("Error, %")
	a1.set_xlabel("Hydraulic diffusivity, m^2/s")
elif file_path_label == 'A':
	a0.set_xlabel("Perpendicular diffusivity, m^2/s")
	a0.set_ylabel("Error, %")
	a1.set_xlabel("Parallel diffusivity, m^2/s")


plot3D = raw_input("3D plot? Y/N:  ").upper()
if plot3D == 'Y':
	fig2 = plt.figure()
	ax = fig2.add_subplot(111, projection='3d')
	ax.set_xlabel('Perpendicular diffusivity (m^2/s)')
	ax.set_ylabel('Parallel diffusivity (m^2/s)')
	ax.set_zlabel('Error (%)')

for constituent in tide_list:
	for case_number in case_list:
		if len(case_list)==1 or choice1=='M':
			min_pts_raw = get_points(output_file, constituent, case_number)
			min_pts = point_to_list(min_pts_raw)
		
		if choice1 == 'O':
# Currently only doing .points file, no .plot file.
			all_pts_raw = get_points(points_file, constituent, case_number)
			all_pts = point_to_list(all_pts_raw)

			x1 = np.array(all_pts[0])
			x2 = np.array(all_pts[1])
			y = np.array(all_pts[2])

		if len(case_list)==1 or choice1=='M':
			min_x1 = np.array(min_pts[0])
			min_x2 = np.array(min_pts[1])
			min_y = np.array(min_pts[2])

		if len(case_list)==1:
			plotlabel = constituent
		elif len(tide_list)==1:
			plotlabel = 'Case '+case_number
		else:
			plotlabel = constituent+','+case_number
		
		if choice1=='M':
			a0.plot(min_x1,min_y, '.', label = plotlabel)
			a1.plot(min_x2,min_y, '.')
					
		elif choice1=='O':
			a0.plot(x1,y, '.', label = plotlabel)
			a1.plot(x2,y, '.')
			if len(case_list)==1:
				a0.plot(min_x1,min_y,'or')
				a1.plot(min_x2,min_y, 'or')

		if plot3D == 'Y':
			if choice1=='M':
				ax.plot(min_x1,min_x2,min_y,'.')		
			elif choice1 == 'O':
				ax.plot(x1,x2,y,'.')
				if len(case_list)==1:
					ax.plot(min_x1,min_x2,min_y,'or')

choice100 = raw_input("Print legend? Y/N:  ").upper()
if choice100 == 'Y':
	a0.legend()

plt.show()

if (len(case_list) > 1) and choice1=='M':
	plt.plot(min_x1, min_x2, '.')
	# This tidbit added Nov 6, 2017; for single constituents and cases only.
	# To get it to plot when it asks for more than one case, enter a case that doesn't exist FIRST,
	# then enter the case you're interested in.
	
	plt.show()
