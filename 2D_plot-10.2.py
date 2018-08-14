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
# v10.2: Created January 1, 2017.
# This is a barebones plotter for the 2D_...-10 stuff. This is to allow me to
#	plot anything outside the usual.
#
#==============================================================================
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np

#==============================================================================
print '\n ======== \n'

# Getting constituents and cases to plot.
tide_list = []
case_raw = []

print "Enter each constituent of interest, one at a time."
print "Enter 'D' when you are done."
while True:
	tide_input = raw_input(">  ").upper()
	if (tide_input == 'D'):
		break
	tide_list.append(tide_input)

print ' '
print "Enter each case of interest, one at a time."
print "Enter 'D' when you are done."
print "Enter 'all' if you would like all cases from the .out file."
while True:
	case_input = raw_input(">  ").upper()
	if (case_input == 'D') or (case_input == 'ALL'):
		break
	case_raw.append(case_input)

file = sys.argv[1]

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


ani_or_iso = raw_input("(I)sotropic, (A)nisotropic, or (O)ther?  ").upper()


#PLOTTING!
fig,(a0,a1) = plt.subplots(1,2, sharey = True)
plot3D = raw_input("3D plot? Y/N:  ").upper()
if plot3D == 'Y':
	fig2 = plt.figure()
	ax = fig2.add_subplot(111, projection='3d')

	if ani_or_iso == 'I':
		ax.set_xlabel('Loading eff.')
		ax.set_ylabel('Diffusivity (m^2/s)')
		ax.set_zlabel('Error (%)')
#	elif ani_or_iso == 'A':
	elif ani_or_iso == 'O':
		ax.set_xlabel('Perpendicular diffusivity (m^2/s)')
		ax.set_ylabel('Parallel diffusivity (m^2/s)')
		ax.set_zlabel('Error (%)')

for constituent in tide_list:
	for case_number in case_list:
# Currently only doing .points file, no .plot file.
		all_pts_raw = get_points(file, constituent, case_number)
		all_pts = point_to_list(all_pts_raw)
		
		x1 = np.array(all_pts[0])
		x2 = np.array(all_pts[1])
		y = np.array(all_pts[2])

		if len(case_list)==1:
			plotlabel = constituent
		elif len(tide_list)==1:
			plotlabel = 'Case '+case_number
		else:
			plotlabel = constituent+','+case_number
		
		a0.plot(x1,y, '.', label = plotlabel)
		a1.plot(x2,y, '.')

		if plot3D == 'Y':
			ax.plot(x1,x2,y,'.', label = plotlabel)
			
	if ani_or_iso == 'I':
#	if ani_or_iso == 'O':
		a0.set_xlabel("Loading efficiency")
		a0.set_ylabel("Error (%)")
		a1.set_xlabel("Hydraulic diffusivity (m^2/s)")
#	elif ani_or_iso == 'A':
	elif ani_or_iso == 'O':
		a0.set_xlabel("Perpendicular diffusivity (m^2/s)")
		a0.set_ylabel("Error (%)")
		a1.set_xlabel("Parallel diffusivity (m^2/s)")

choice100 = raw_input("Print legend? Y/N:  ").upper()
if choice100 == 'Y':
	a0.legend()
	if plot3D == 'Y':
		ax.legend()

plt.show()


