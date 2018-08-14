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
# On the prompt line: python 2D_plot-10.0.py <file_path> <file_title (exclude extension)>
#
#==============================================================================
# v10.1: Created October 13, 2016.
#
# This script includes how to deal with the points file included in 2D_anisotropic-10.1.py
#
#==============================================================================
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import numpy as np

#==============================================================================
print '\n ======== \n'

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

tide = raw_input("What is the constituent?  ")
case = raw_input("What is the case number?  ")
file_path = sys.argv[1]
file_title = sys.argv[2]

if file_path.upper() == 'A':
	file_path = "/Users/admin/desktop/tidal analysis/archives/2d_anisotropic-10/"
elif file_path.upper() == 'I':
	file_path = "/Users/admin/desktop/tidal analysis/archives/2d_isotropic-10/"

choice1 = raw_input("p(L)ot or p(O)ints?  ").upper()

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
			if ("Source".upper() in LINE):
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
min_points = []

min_pts_raw = get_points(output_file, tide, case)

if choice1 == 'L':
	all_pts_raw = get_points(plot_file, tide, case)
	#print all_pts_raw
elif choice1 == 'O':
	all_pts_raw = get_points(points_file, tide, case)
	#print all_pts_raw

min_pts = point_to_list(min_pts_raw)
all_pts = point_to_list(all_pts_raw)
#print all_pts

x1 = np.array(all_pts[0])
x2 = np.array(all_pts[1])
y = np.array(all_pts[2])

min_x1 = np.array(min_pts[0])
min_x2 = np.array(min_pts[1])
min_y = np.array(min_pts[2])

fig,(a0,a1) = plt.subplots(1,2, sharey = True)
a0.plot(x1,y, '.')
a0.plot(min_x1,min_y,'or')
a1.plot(x2,y, '.')
a1.plot(min_x2,min_y, 'or')

plot3D = raw_input("3D plot? Y/N:  ").upper()
if plot3D == 'Y':
	fig2 = plt.figure()
	ax = fig2.add_subplot(111, projection='3d')
	ax.plot(x1,x2,y,'.')
	ax.plot(min_x1,min_x2,min_y,'or')

plt.show()