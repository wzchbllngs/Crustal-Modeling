# 1. Isotropic or anisotropic
# 2. get number of sources
# 3. get details of each source
# 4. calculate "tide" signal/amplitudes

# This script generates input files of synthetic cases to test 2D_isotropic-10
#	and 2D_anisotropic-10 models. It generates them based on the ideal
#	isotropic and anisotropic model(s) I've been looking/dealing with.
#
# Version 10.1: April 21, 2018
#	This version includes a term for the earth tide in the ratio_phi() function.
# ===========================================================================80
import numpy as np

# ===========================================================================80
# Givens:

tidal_periods = {'M2': 44714.
				}
tides = tidal_periods.keys()

# Location for .input file:
default_out = "/Users/admin/desktop/Tidal Analysis/Archives/2D_test-10/"

# Location for tidal data file:
default_data = "/Users/admin/desktop/Tidal Analysis/1234ABCD/Tidal Data/"
suffix = ".input"

sigma = 8. # Seafloor amplitude
phi_sea = 0. # seafloor phase

earth_amplitude = 0.01
phi_earth = -np.pi/3. #(basically -60 degrees is out of the range of the loading signal)

# Parameters for .input file:
hole = "1234ABCD"
gammastep = .01
etastep = 5. #m^2/s
etacap = 3000. #m^2/s

# ===========================================================================80
# User inputs:

print ' '
print " ~~~~~~~~~~~~~~~~~~~~"
print ' '

print default_out
#path_choice = raw_input("(D)efault output path, or (C)ustom?  ").lower()
path_choice = 'd'
while True:
	if path_choice=='d' or path_choice=='c':
		break
	path_choice = raw_input("Please make a valid choice:  ").lower()

if path_choice=='d':
	output_path = default_out
if path_choice=='c':
	output_path = raw_input("Please enter output path now:  ")

if output_path[-1] != '/':
	output_path += '/'

filename = raw_input(".input file name?  ")

# -------------------------------------

print " "
print tides
#print "Which would you like to exclude, if any?"
#print "(Enter 'D' when done)"
#tide_choice = raw_input(">  ").upper()
tide_choice = 'D'
while True:
	ctr = 0
	if tide_choice == 'D':
		break
	if tide_choice not in tides:
		print "Please make a valid selection."
	if len(tides)==0:
		print "You've removed all constituents!"
		exit()
	for tide in tides:
		if tide_choice==tide.upper():
			tides.pop(ctr)
		ctr += 1
	tide_choice = raw_input(">  ").upper()

print tides

print " "
print " ~~~~ Parameters:"
print "Which model?"
print "1. Isotropic"
print "2. Anisotropic"

#model = raw_input(">  ")
model = '2'
while True:
	if model=='1' or model=='2':
		break
	model = raw_input("Please make a valid selection:  ")
if model=='1':
	output_path += "Isotropic/"
elif model=='2':
	output_path += "Anisotropic/"

# ~~~~ Taking care of the sources:
source_dict = {}
print ' '
no_sources = int(raw_input("How many sources?  "))
while True:
	if no_sources > 0:
		break
	no_sources = int(raw_input("Please enter a valid number:  "))

ctr = 1
while ctr <= no_sources:
	print "Source "+str(ctr)+":"
	distance = float(raw_input("Distance from site in meters:  "))
#	gamma = float(raw_input("Loading efficiency at source (default=1):  "))
	gamma = 1.

	if model=='1':
		source_dict[ctr] = (distance, gamma)

	if model=='2':	# Anisotropic
		angle = float(raw_input("Angle (deg) from perpendicular:  "))
		source_dict[ctr] = (distance, angle, gamma)
	print ' '
	ctr += 1

# ~~~~ Aquifer properties:
print " ~~~~ Aquifer properties:"

#gamma_site = float(raw_input("Loading efficiency at site:  "))
gamma_site = 0.6

if model=='1':		# Isotropic
	eta = float(raw_input("Hydraulic diffusivity (m^2/s):  "))
#	eta = 1000.
elif model=='2':	# Anisotropic
	print "Hydraulic diffusivity in m^2/s..."
#	eta_x = float(raw_input("...along perpendicular:  "))
#	eta_y = float(raw_input("...along parallel:  "))
	eta_x = 1000.
	eta_y = 500.
	eta = (eta_x, eta_y)

# ===========================================================================80
# Functions:

# Outputs the phase lag and ratio of the tide signal.
def ratio_phi(gamma_site,eta,period, source_dict, model):
	# model: 1 = isotropic, 2 = anisotropic
	
	diff_cos = 0
	diff_sin = 0
	
	earth_amplitude = 0.01
	phi_earth = -np.pi/180. * -60. #(basically negative degrees is out of the range of the loading signal)

	
	for source in source_dict:
		if model==str('1'):
			theta = 0
			gamma_src = source_dict[source][1]
			eta_x = eta
			eta_y = eta
		
		if model==str('2'):
			theta = source_dict[source][1]*np.pi/180.
			gamma_src = source_dict[source][2]
			eta_x = eta[0]
			eta_y = eta[1]
		
		z = source_dict[source][0]
		eta_eff = eta_x*(np.cos(theta))**2 + eta_y*(np.sin(theta))**2
		Lambda = np.pi*z / np.sqrt(np.pi*eta_eff*period)
		
		
		
		cos = (gamma_src - gamma_site) * np.exp(-Lambda) * np.cos(Lambda)
		sin = (gamma_src - gamma_site) * np.exp(-Lambda) * np.sin(Lambda)
		diff_cos += cos
		diff_sin += sin
	
	cos_term = diff_cos + gamma_site + earth_amplitude*np.cos(phi_earth)
	sin_term = diff_sin + earth_amplitude*np.sin(phi_earth)
	
	ratio = np.sqrt(sin_term**2 + cos_term**2)
	phase_lag = np.arctan(sin_term/cos_term)
	
	return(ratio, phase_lag)

# ===========================================================================80
# Creating files:

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Creating data files in '.../Tidal Analysis/1234ABCD/' ~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

borehole_file = open(default_data+filename+' borehole.txt', 'w')
seafloor_file = open(default_data+filename+' seafloor.txt', 'w')

# ~~~~ Writing down aquifer parameters:	
borehole_file.write("# Aquifer parameters:\n")
borehole_file.write("# ===================\n")
seafloor_file.write("# Aquifer parameters:\n")
seafloor_file.write("# ===================\n")
if model=='1': # Isotropic
	borehole_file.write("#  Isotropic\n")
	seafloor_file.write("#  Isotropic\n")
elif model=='2': # Anisotropic
	borehole_file.write("#  Anisotropic\n")
	seafloor_file.write("#  Anisotropic\n")

borehole_file.write("#    Loading efficiency @ site: {}\n".format(gamma_site))
seafloor_file.write("#    Loading efficiency @ site: {}\n".format(gamma_site))
if model=='1':
	borehole_file.write("#    Hydraulic diffusivity (m^2/s): {}\n".format(eta))
	seafloor_file.write("#    Hydraulic diffusivity (m^2/s): {}\n".format(eta))
elif model=='2':
	borehole_file.write("#    Hydraulic diffusivity (m^2/s):\n")
	borehole_file.write("#      Perpendicular to axis: {}\n".format(eta[0]))
	borehole_file.write("#      Parallel to axis: {}\n".format(eta[1]))
	seafloor_file.write("#    Hydraulic diffusivity (m^2/s):\n")
	seafloor_file.write("#      Perpendicular to axis: {}\n".format(eta[0]))
	seafloor_file.write("#       Parallel to axis: {}\n".format(eta[1]))
	
borehole_file.write("#\n")
seafloor_file.write("#\n")
borehole_file.write("# Sources\n")
borehole_file.write("# =======\n")
seafloor_file.write("# Sources\n")
seafloor_file.write("# =======\n")
if model=='1':
	borehole_file.write("#  Source # = (distance [m], loading eff. @ source)\n")
	seafloor_file.write("#  Source # = (distance [m], loading eff. @ source)\n")
elif model=='2':
	borehole_file.write("#  Source # = (distance [m], angle [deg], loading eff. @ source)\n")
	seafloor_file.write("#  Source # = (distance [m], angle [deg], loading eff. @ source)\n")

for source in source_dict:
	borehole_file.write("#    Source {} = {}\n".format(source,source_dict[source]))
	seafloor_file.write("#    Source {} = {}\n".format(source,source_dict[source]))

# ~~~~ Writing down the tidal data:
borehole_file.write("#\n")
seafloor_file.write("#\n")
for constituent in tides:
	period = tidal_periods[constituent]
	raw_data = ratio_phi(gamma_site,eta,period , source_dict,model)
	speed = 360./(period/3600.) # degrees/hour
	
#	earth_signal = earth_amplitude
	p = sigma * raw_data[0] # borehole pressure
	phi_bore = raw_data[1] * 180./np.pi

	borehole_file.write('{:>10} {:>12.8f} {:>12} {:12.4f}\n'.format(constituent,speed,p,phi_bore))
	seafloor_file.write('{:>10} {:>12.8f} {:>12} {:12.4f}\n'.format(constituent,speed,sigma,phi_sea))
	
borehole_file.close()
seafloor_file.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Creating .input file for 2D_isotropic-10 or 2D_anisotropic-10 ~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

input_file = open(output_path+filename+'.input','w')

# ~~~~ Writing down aquifer parameters
input_file.write("# Aquifer parameters:\n")
input_file.write("# ===================\n")

if model=='1': # Isotropic
	input_file.write("#  Isotropic\n")
elif model=='2': # Anisotropic
	input_file.write("#  Anisotropic\n")

input_file.write("#    Loading efficiency @ site: {}\n".format(gamma_site))
if model=='1':
	input_file.write("#    Hydraulic diffusivity (m^2/s): {}\n".format(eta))
elif model=='2':
	input_file.write("#    Hydraulic diffusivity (m^2/s):\n")
	input_file.write("#      Perpendicular to axis: {}\n".format(eta[0]))
	input_file.write("#      Parallel to axis: {}\n".format(eta[1]))
	
input_file.write("#\n")
input_file.write("# Sources\n")
input_file.write("# =======\n")
if model=='1':
	input_file.write("#  Source # = (distance [m], loading eff. @ source)\n")
elif model=='2':
	input_file.write("#  Source # = (distance [m], angle [deg], loading eff. @ source)\n")

for source in source_dict:
	input_file.write("#    Source {} = {}\n".format(source,source_dict[source]))

input_file.write("#\n")
input_file.write("hole = {};\n".format(hole))
input_file.write("dataset = {};\n".format(filename))

input_file.write("datapath = /Users/admin/desktop/Tidal Analysis/1234ABCD/Tidal Data/;\n")
if model=='1':
	input_file.write("outputpath = /Users/admin/desktop/Tidal Analysis/Archives/" + \
"2D_test-10/Isotropic/;\n")
elif model=='2':
	input_file.write("outputpath = /Users/admin/desktop/Tidal Analysis/Archives/" + \
"2D_test-10/Anisotropic/;\n")

for constituent in tides:
	if model=='1':
		input_file.write("tide = {};\n".format(constituent))
	elif model=='2':
		input_file.write("tide = {} {};\n".format(constituent, gamma_site))
	
for source in source_dict:
	input_file.write("Source "+str(source)+" =")
	for item in source_dict[source]:
		input_file.write(' {}'.format(item))
	input_file.write(';\n')

if model=='1':
	input_file.write("gammastep = {};\n".format(gammastep))
input_file.write("etastep = {};\n".format(etastep))
input_file.write("etacap = {};\n".format(etacap))
	


