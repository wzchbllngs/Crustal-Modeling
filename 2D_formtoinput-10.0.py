# This script takes a .form file (title of the form <site>.form) and creates
#	a specified number of .input files from it.
#
#==============================================================================
import os
import itertools

path0 = "/Users/admin/desktop/Tidal analysis/Archives/"
pathIsotropic = path0 + '2D_isotropic-10/'
pathAnisotropic = path0 + '2D_anisotropic-10/'

alphabet = ('a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',
'p','q','r','s','t','u','v','w','x','y','z')
#==============================================================================
# FUNCTIONS:
	

#==============================================================================
print '\n========\n'

# Gathering basic data about what the user wants.

hole = raw_input("Site?  ").upper()
wanted_file = hole+'.form'

print "From the default paths..."
choice = raw_input("(I)sotropic, (A)nisotropic, or (C)ustom path?  ").upper()
while True:
	if choice=='I' or choice=='A' or choice=='C':
		break
	choice = raw_input("Please enter a valid selection:  ").upper()

if choice=='I':
	path = pathIsotropic
elif choice=='A':
	path = pathAnisotropic
elif choice=='C':
	print "Please enter the path your .form file is located in."
	path = raw_input(">  ")
	if path[-1] != '/':
		path += '/'

data_list = []
print "Years/datasets to include? Enter 'D' when done."
while True:
	set = raw_input(">  ")
	if set.upper() == 'D':
		break
	data_list.append(set)

# Getting the numbers and labels of each source from the .form file.

## Checking to see if the file is in the given path:
directory = os.listdir(path)
if wanted_file not in directory:
	print wanted_file + " not found in"
	print path
	exit()

## Getting sources and storing the distances of each:
lines = []

sources = []
source_dict = {}
file_form = open(path+wanted_file, 'r')
for line in file_form:
	if line=='\n':
		continue
	if line[0] != '#':
		print line
		if line[:6].lower()=='source':
			linesplit = line.split(' ')
			source = linesplit[1]
			end = linesplit[-1]
			endsplit = end.split(';')
			distance = float(endsplit[0])
			if source not in sources:
				sources.append(source)
				source_dict[source] = []
			if distance not in source_dict[source]:
				source_dict[source].append(distance)
		else:
			lines.append(line)

file_form.close()
			
for source in sources:
	source_dict[source].sort()
print sources
print source_dict
print lines

# Creating the .input files.
## Gives choice over whether to make all combos of sources, or
##	to just go in order from closest to farthest.

print "Would you like to..."
print "1. create (A)ll combos of sources, or"
print "2. or (G)roup from closest to farthest?"
choice2 = raw_input(">  ").upper()
while True:
	if choice2=='1' or choice2=='A':
		choice2 = 'A'
		break
	if choice2=='2' or choice2=='G':
		choice2 = 'G'
		break
	choice2 = raw_input("Please make a valid selection.  ").upper()

if choice2=='G':
	print "Do you want the order based on:"
	print "1. Closest possible distances,"
	print "2. Farthest possible distances, or"
	print "3. The average of the two distances?"
	choice3 = raw_input(">  ")
		
	while True:
		if choice3=='1' or choice3=='2' or choice3=='3':
			break
		choice3 = raw_input("Please make a valid selection:  ")



for dataset in data_list:
	no_sources = 1

## Creates files where there is only one source per file.
### ctr keeps track of which file for each total of sources included per .input
###		file--makes the labels a,b,c,d,e...for each file.
	if no_sources==1:
		ctr = 0
		for source in sources:
			title = hole+','+dataset+'--'+str(no_sources)+alphabet[ctr]+'.input'
			file_input = open(path+title, 'w')

			file_input.write("dataset = "+dataset+";\n")
			for line in lines:
				file_input.write(line)
			file_input.write('\n\n')
			for distance in source_dict[source]:
				file_input.write("Source "+source+" = "+str(distance)+";\n")
			file_input.close()
			ctr += 1
	no_sources += 1

### Creating .input files with more than one source for file.
	if choice2=='A':
		while no_sources <= len(sources):
			source_combos = list(itertools.combinations(sources,no_sources))
			ctr = 0
			for combo in source_combos:
				title = hole+','+dataset+'--'+str(no_sources)+alphabet[ctr]+'.input'
				file_input = open(path+title, 'w')
				file_input.write("dataset = "+dataset+";\n")
				for line in lines:
					file_input.write(line)
				file_input.write('\n\n')
				for source in combo:
					for distance in source_dict[source]:
						file_input.write("Source "+source+" = "+str(distance)+";\n")
					file_input.write('\n')
				file_input.close()
				ctr += 1
			no_sources += 1

	elif choice2=='G':
		ordered_distances = []
		ordered_sources = []
		if choice3=='1':
			for source in sources:
				ordered_distances.append(source_dict[source][0])
			ordered_distances.sort()
			if len(ordered_distances) != len(sources):
				print "ordered_distances != sources!"
				exit()
			ii = 0
			while ii < len(sources):
				distance_to_compare = ordered_distances[ii]
				for source in sources:
					if source_dict[source][0]==distance_to_compare:
						if source not in ordered_sources:
							ordered_sources.append(source)
				ii += 1
		
		elif choice3=='2':
			for source in sources:
				ordered_distances.append(source_dict[source][-1])
			ordered_distances.sort()
			if len(ordered_distances) != len(sources):
				print "ordered_distances != sources!"
				exit()
			ii = 0
			while ii < len(sources):
				distance_to_compare = ordered_distances[ii]
				
				for source in sources:
					if source_dict[source][-1]==distance_to_compare:
						if source not in ordered_sources:
							ordered_sources.append(source)
				ii += 1

		elif choice3=='3':
			source_avgs = {}
			for source in sources:
				if len(source_dict[source])==1:
					source_avgs[source] = source_dict[source][0]
				else:
					min = source_dict[source][0]
					max = source_dict[source][-1]
					avg = (min+max)/2.
					source_avgs[source] = avg
			
			for source in sources:
				ordered_distances.append(source_avgs[source])
			
			ordered_distances.sort()
			if len(ordered_distances) != len(sources):
				print "ordered_distances != sources!"
				exit()
			ii = 0
			while ii < len(sources):
				distance_to_compare = ordered_distances[ii]
				
				for source in sources:
					if source_avgs[source]==distance_to_compare:
						if source not in ordered_sources:
							ordered_sources.append(source)
				ii += 1

		if len(ordered_sources) != len(sources):
			print "len(ordered_sources) != len(sources)!"
			exit()		
#		print ordered_sources
#		print ordered_distances
	
		no_sources = 2
		source_combo = [ordered_sources[0]]
		while no_sources <= len(sources):
			source_combo.append(ordered_sources[no_sources-1])
			title = hole+','+dataset+'--'+str(no_sources)+'.input'
			
			file_input = open(path+title, 'w')
			file_input.write("dataset = "+dataset+";\n")
			for line in lines:
				file_input.write(line)
			file_input.write('\n\n')
			
			for source in source_combo:
				for distance in source_dict[source]:
					file_input.write("Source "+source+" = "+str(distance)+";\n")
				file_input.write('\n')
			
			file_input.close()
			no_sources += 1

