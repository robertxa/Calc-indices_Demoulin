######!/usr/bin/env python
# -*- coding: utf-8 -*-
###########################################
#                                         #
#       Calcul Indices Demoulin           #
#         by Thomas Croissant             #
#          and Xavier Robert              #
#           April-mai 2013                #
#                                         #
###########################################


# This script is dedicated to extract Demoulin indices by using the equations
# from Demoulin 2011 and 2012.
# You need to run this script after the ArcGIS script also provided,
# that generate .dxf files, base of this calculation
# This srcipt generate :
#	- for each river, a .txt file with all the calculations stored in
#	  rivers folder, and hypsometric curves stored in Graphs folder
#	- a general .txt file summarizing the cacul of the indices stored in the working folder
#	- a graph representing the regression line for the system studied
#
# You also need to write a file "River-caract.txt" that has this structure :
#	- Line 1 = header (river_nb	L_watershed)
#	- Line 2-->n+1 = i	length of the watershed for the river i (i:1-->n)
#	- column should be delimited with a tab (\t)
#
# To keep a clean workspace, I built this set of scripts on a strutured project.
# You need to follow this structure for the database/.project:
#	|_RASTER
#	|	|_raster_fnme+"_corr"
#	|	|_raster_fnme+"_fdir"
#	|	|_raster_fnme+"_facc"
#	|	|_raster_fnme+"_alti"
#	|	|_raster_fnme+"_cost"
#	|_RPROFILS
#	|	|_Ri
#	|	|	|_TMP
#	|	|	|	|_tem_files_produced_by_this_script
#	|	|	|_SHP_FILES
#	|	|	|	|_exu_Ri
#	|	|	|	|_src_Ri
#	|_script_ArcGIS.py
#	|_Calcul-indices_Demoulin.py
#	|_River-caract.txt
#
# Read the comments in the script for more informations.




# Import packages needed. Be careful with versions !
import datetime
from dbfpy import dbf                            # to read dbf tables http://dbfpy.sourceforge.net
import math                                      # For calculations
import numpy as np                               # need version 1.7 or higher
import scipy as sp                               # need version 0.12 or higher
from scipy import integrate, stats               # to calculate integrals and regressions
import matplotlib.pyplot as plt                  # module to plot figures
from pylab import savefig                        # Module to manage figures
from os import path, access, R_OK, mkdir         # W_OK for write permission.
import sys                                       # Module to use system commands
import copy										 # To copy variables


###############################################################################
###############################################################################
####### SET your ENVIRONNEMENT  #################
# Declaration of variables
# You can edit them for your workspace settings
# Pixel espacement in the dem (in m)
dnx = 100
dny = 100
factor = 1000000
# Files containing the caracteristic of the rivers you want to run:
#     first line =  header [river_nb,L_watershed]
#     then, one line per river, with 2 columns [river_nb,L_watershed]
#     tab separation between the columns
Rivercaract = "River-caract.txt"
# Folder where will be stored the graphs
graphs_path = "Graphs"                             
# Folder where are stored the rivers data
rprofiles = "RProfils/"
######## END of modifications ###################                            
###############################################################################
###############################################################################

print "Script to extract data to apply Demoulin approach"
print "on rivers profiles"
print "X. Robert, Grenoble, 06/2013"
print "  "
print "  "


# Check the versions of numpy and scipy, and return error message if not compatible and quit
if np.__version__ < '1.7.1'and sp.__version__ <'0.12':
	 sys.exit('ERROR : numpy version should be 1.7 or higher and scipy version \
		  	   should be 0.12 or higher ; Upgrade them')  
if np.__version__ < '1.7.1'and sp.__version__ >='0.12':
	 sys.exit('ERROR : numpy version should be 1.7 ; Upgrade it')
if np.__version__ >= '1.7.1'and sp.__version__ <'0.12':
	 sys.exit('ERROR : scipy version should be 0.12 or higher ; Upgrade it') 	  	 

#generate arrays to store datas from files and calculs  
  # calcul_indices.xls
px_cum = []
Lwatershed = []
riversnb = []
Efact = []

# Verify if the folder for Graphs exists
if path.exists(graphs_path) == False:
	# If it do not exists, create it
	mkdir(graphs_path)
  
# read the file River-charact.txt
# Should contain Rivers numbers with length of the watershed. See example provided
print ('Reading {FileNa}'.format(FileNa=str(Rivercaract)))
if path.isfile(Rivercaract) and access(Rivercaract, R_OK):            #Check if the file exist
	river = open(Rivercaract, "r")
	nd = 0                                                            #nd = number of rivers to read
	# Read Header, remove end of line, and extract the fields separated by a Tab 
	entete = river.readline().rstrip('\n\r').split('\t')
	# Determine index of the different fields
	rivernbidx = entete.index("river_nb")
	Lwatershedidx = entete.index("L_watershed")
	for line in river:                                                              
		datar = line.rstrip('\n\r').split('\t')                       #read the line, separated with tab (\t)
		riversnb.append(datar[rivernbidx])
		Lwatershed.append(float(datar[Lwatershedidx]))
		nd += 1	                                                      #read the total number of lines/rivers
	print "Number of rivers :", nd
	print "   "
else:                                                                 #If file does not exist, print error message and quit
	print('ERROR : File {FileNa} does not exist'.format(FileNa=str(Rivercaract))) 
	sys.exit('File {FileNa} should contain an hearder (first line), and then the number \
			  of the river and the length of the watershed'.format(FileNa=str(Rivercaract))) 

# Generate numpy arrays
Ib=np.empty(nd)
Ir=np.empty(nd)
lnA=np.empty(nd)
IrIb=np.empty(nd)
Efact=np.empty(nd)
R1w=np.empty(nd)

f0w = open("Summary.txt", "w") # writting
line = "Summary of Calcul-indices-Demoulin.py script \n"
f0w.write(line)
line = "by Xavier Robert, Grenoble, 06/2013 \n"
f0w.write(line)
line = " \n"
f0w.write(line)
f0w.write(line)
line = "Number of rivers :" + str(nd) +"\n"
line = " \n"
f0w.write(line)

# Do a loop on the number of Rivers
for k in range(nd):
	# Init graph
	plt.clf()
	print "Working on river number ",k+1
	line = "River number" + str(k+1) + "\n"
	f0w.write(line)
	#generate arrays to store datas from files and calculs
	# Rk.xls
	pxvalue = []            # river pixel value data from R1.dbf
	pxcount = []            # river pixel count data from R1.dbf
	pxmin = []              # min pixel value data from R1.dbf
	pxmax = []              # max pixel value data from R1.dbf
	alti_bass = []          # (watershed) from R1demb.dbf
	px_bass = []            # (watershed) from R1demb.dbf
	px_bass_cumul = []
	alti_network = []       # from R1demr.dbf
	px_network = []         # from R1demr.dbf

	# Read data from the input files --> to vectors
	fnme = "R"+str(k+1)
	data_file=fnme+".dbf"                             #add the extention for the txt file
	data_file_path=rprofiles+fnme+"/"+data_file       #change name in path
	print data_file_path
	print "Reading dbf files",fnme+"/"+fnme+".dbf"
	if path.exists(data_file_path) and path.isfile(data_file_path) and access(data_file_path, R_OK):       #Check if the file exist
		f1i = dbf.Dbf(rprofiles+fnme+"/"+fnme+".dbf")                               # reading dbf
		rdr1 = 0
		for line in f1i:
			# table structure :
			# VALUE; COUNT; AREA; MIN; MAX; RANGE; MEAN; STD; SUM; VARIETY; MAJORITY; MINORITY; MEDIAN		
			pxvalue.append(line[0])		
			pxcount.append(line[1])
			pxmin.append(line[3])
			pxmax.append(line[4])
			rdr1 += 1
	else:                                                                  #If file does not exist, print message and quit
		sys.exit('ERROR : File {FileNa} does not exist'.format(FileNa=str(fnme+"/"+fnme+".dbf"))) 	
	print "Reading dbf files",fnme+"/"+fnme+"demb.dbf"	
	data_file=fnme+"demb.dbf"                             #add the extention for the txt file
	data_file_path=rprofiles+fnme+"/"+data_file           #change name in path
	if path.exists(data_file_path) and path.isfile(data_file_path) and access(data_file_path, R_OK):       #Check if the file exist
		f2i = dbf.Dbf(rprofiles+fnme+"/"+fnme+"demb.dbf")                  # reading dbf
		for line in f2i:                                                   # Skip the first one		
			alti_bass.append(float(line[0]))
			px_bass.append(float(line[1]))
		rdr2 = len(px_bass)
	else:                                                                  #If file does not exist, print message and quit
		sys.exit('ERROR : File {FileNa} does not exist'.format(FileNa=str(fnme+"/"+fnme+"demb.dbf"))) 
	print "Reading dbf files",fnme+"/"+fnme+"demb.dbf"
	data_file=fnme+"demr.dbf"                            #add the extention for the txt file
	data_file_path=rprofiles+fnme+"/"+data_file          #change name in path
	if path.exists(data_file_path) and path.isfile(data_file_path) and access(data_file_path, R_OK):       #Check if the file exist
		f3i = dbf.Dbf(rprofiles+fnme+"/"+fnme+"demr.dbf")                  # reading dbf
		for line in f3i:                                                   # Skip the first one
			alti_network.append(float(line[0]))
			px_network.append(float(line[1]))
		rdr3 = len(px_network)
	else:                                                                  #If file does not exist, print message
		sys.exit('ERROR : File {FileNa} does not exist'.format(FileNa=str(fnme+"/"+fnme+"demr.dbf"))) 
		
		# close the 3 files
	f1i.close()
	f2i.close()
	f3i.close()
		
	print "Calcul table"
	# Calculate first table
	# extract min and max values
	maxa = float(pxvalue[rdr1-1])
	ming = float(alti_bass[0])
	maxg = float(alti_bass[rdr2-1])

	ln=max(rdr1,rdr2,rdr3)                                              # Calcul max length of the vectors
	# Generate numpy arrays
	val_norm = np.empty(len(pxmin))
	min_norm = np.empty(len(pxmin))
	for i in range(rdr1):
		val_norm[i] = pxvalue[i]/maxa
		min_norm[i] = (pxmin[i]-ming)/(maxg-ming)
	
	# Generate numpy arrays
	px_bass_cumul = np.empty(len(px_bass))
	alti_bass_norm = np.empty(len(px_bass))
	px_bass_norm = np.empty(len(px_bass))	
	px_bass_cumul[rdr2-1] = 1.0
	for i in range(rdr2):
		if i == 0:
			px_bass_cumul[rdr2-1] = 1.0
		else:
			px_bass_cumul[rdr2-i-1] = px_bass_cumul[rdr2-i] + px_bass[rdr2-i-1]
		alti_bass_norm[i] = (alti_bass[i]-ming)/(maxg-ming)
	maxI = float(px_bass_cumul[0])
	for i in range(rdr2):  		
		px_bass_norm[i] = px_bass_cumul[i]/maxI
  	
  	# Generate numpy arrays	
	px_network_cumul = np.empty(len(px_network))
	alt_network_norm = np.empty(len(alti_network))
	px_network_norm = np.empty(len(alti_network))
	for i in range(rdr3):
		if i == 0:
			px_network_cumul[rdr3-1] = 1.0
		else:
			px_network_cumul[rdr3-i-1] = px_network_cumul[rdr3-i] + px_network[rdr3-i-1]
		alt_network_norm[i] = (alti_network[i]-ming)/(maxg-ming)		
	maxn = float(px_network_cumul[0])
	for i in range(rdr3):
		px_network_norm[i] = px_network_cumul[i]/maxn
	
	# write the new table
	f1w = open(rprofiles+fnme+"/R"+str(k+1)+"-calc.txt", "w") # writting
	# write header
	header = "VALUE \t COUNT \t MIN \t MAX \t Val_norm \t min_norm \t alti_bass \t \
	          px_bass \t px_bass_cumul \t alti_bass_norm \t px_bass_norm \t alt_network \
	          \t px_network \t px_network_cum \t alti_network_norm \t px_network_norm \n"
	f1w.write(header)
	
	# give the same length to each list regardings the longest :	
	if rdr1<ln:
		pxvalue[rdr1:] = ["NaN"]*(ln-rdr1)
		pxcount[rdr1:] = ["NaN"]*(ln-rdr1)
		pxmin[rdr1:] = ["NaN"]*(ln-rdr1)
		pxmax[rdr1:] = ["NaN"]*(ln-rdr1)
		# sequence to copy an rdr1 array to a longest array that we finish to fill with NaN Values
		# Method 1)
		#	d=copy.copy(c) #: comes from copy module
		#	d=np.resize(d,[ln]) #: Need numpy
		# Method 2)
		#	temp = new array with a length (ln-rdri)
		#	fill temp with "NaN" value
		#	append temp to d
		temp = np.empty(ln-rdr1)
		temp[0:] = ["NaN"]*(ln-rdr1)
		val_norm_write = np.empty(ln)
		val_norm_write = np.append(val_norm,temp)
		min_norm_write = copy.copy(min_norm)
		min_norm_write = np.resize(min_norm_write,[ln])
		min_norm_write[rdr1:] = ["NaN"]*(ln-rdr1)
	else:
		#if ln = rdr1, copy arrays to write it in the txt file
		val_norm_write = val_norm
		min_norm_write = min_norm	
	
	if rdr2<ln:
		alti_bass[rdr2:] = ["NaN"]*(ln-rdr2)
		px_bass[rdr2:] = ["NaN"]*(ln-rdr2)
		px_bass_cumul[rdr2:] = ["NaN"]*(ln-rdr2)	
		alti_bass_norm_write = copy.copy(alti_bass_norm)
		alti_bass_norm_write = np.resize(alti_bass_norm_write,[ln])
		alti_bass_norm_write[rdr2:] = ["NaN"]*(ln-rdr2)
		px_bass_norm_write = copy.copy(px_bass_norm)
		px_bass_norm_write = np.resize(px_bass_norm_write,[ln])
		px_bass_norm_write[rdr2:] = ["NaN"]*(ln-rdr2)
	else:
		alti_bass_norm_write = alti_bass_norm
		px_bass_norm_write = px_bass_norm

	if rdr3<ln:
		alti_network[rdr3:] = ["NaN"]*(ln-rdr3)
		px_network[rdr3:] = ["NaN"]*(ln-rdr3)		
		temp = np.empty(ln-rdr3)
		temp[0:] = ["NaN"]*(ln-rdr3)
		px_network_cumul_write = np.empty(ln)
		px_network_cumul_write = np.append(px_network_cumul,temp)
		alt_network_norm_write=np.empty(ln)
		alt_network_norm_write = np.append(alt_network_norm,temp)
		px_network_norm_write = copy.copy(px_network_norm)
		px_network_norm_write = np.resize(px_network_norm_write,[ln])
		px_network_norm_write[rdr3:] = ["NaN"]*(ln-rdr3)
	else:
		alt_network_norm_write = alt_network_norm
		px_network_norm_write = px_network_norm

	# write the data in the file under the header			
	for i in range(ln):
		# build the line for each iteraction
		line = str(pxvalue[i]) + "\t" + str(pxcount[i]) + "\t" + str(pxmin[i]) + \
		        "\t" + str(pxmax[i]) + "\t" + str(val_norm_write[i]) + "\t" + \
		        str(min_norm_write[i]) + "\t" + str(alti_bass[i]) + "\t" + \
		        str(px_bass[i]) + "\t" + str(px_bass_cumul[i]) + "\t" + \
		        str(alti_bass_norm_write[i]) + "\t" + str(px_bass_norm_write[i]) +  \
		        "\t" + str(alti_network[i]) + "\t" + str(px_network[i]) + "\t" + \
		        str(px_network_cumul_write[i]) + "\t" + str(alt_network_norm_write[i]) \
		        + "\t" + str(px_network_norm_write[i]) + "\n"
		# write the line
		f1w.write(line)	
	# Close the file			
	f1w.close()

	# Build last table : one line for each river	
	# write in the file Calc-indices-fastSR.txt
	if k == 0:
		f2w = open("Calc-indices-fastSR.txt", "w") # writing
		# write the header when openning the file with the first river computed
		header = "px_cumul \t  Ln(A)_(km2) \t Ib \t Ir \t Ir/Ib \t L_(km) \t E \t R1w \t River_nb \n"
		f2w.write(header)
	print "Calculate hypsometric integrals"


	# Calculate Hypsometric integrals ; need scipy version 0.12 or higher
	hb = sp.integrate.simps(px_bass_norm,alti_bass_norm)
	hn = sp.integrate.simps(px_network_norm,alt_network_norm)
	hr = sp.integrate.simps(min_norm,val_norm)    

	Ib[k]=hb-hn
	Ir[k]=hn-hr
	print "river ",k+1, " hb :",hb,"; hn :",hn,"; hr :",hr, "Ir :", Ir[k], "Ib :", Ib[k]
	line = "river " + str(k+1) + " hb :" + str(hb) + "; hn :" + str(hn) + "; hr :" + \
	       str(hr) + "Ir :" + str(Ir[k]) + "Ib :" + str(Ib[k]) + " \n"
	f0w.write(line)
	line = "  \n"
	f0w.write(line)
	
	# Calcul the Efact and R1w paramters for each river
	lnA[k]=math.log(maxI * dnx * dny / factor)
	IrIb[k]=Ir[k] / Ib[k]
	Efact[k] = 4.0 * math.exp(lnA[k]) / (math.pi * Lwatershed[k]**2)	
	R1w[k] = IrIb[k] / math.sqrt(Efact[k])
	
	# Built the line to be print in the general file
	line = str(maxI) + "\t" + str(lnA[k]) + "\t" + str(Ib[k]) + "\t" + str(Ir[k]) + \
	       "\t" + str(IrIb[k]) + "\t" + str(Lwatershed[k]) + "\t" + str(Efact[k]) + \
	       "\t" + str(R1w[k]) + "\t R" + str(k+1) +"\n"
	# write the line in the general file
	f2w.write(line)
	
	# plot Hysometric curves
	#w, h = figaspect(1.5)
  	#pyplot.figure(figsize=(10,17))
	Title_string="Hypsometric_curves_R"+str(k+1)	
	plt.plot(alti_bass_norm,px_bass_norm,alt_network_norm,px_network_norm,val_norm,min_norm,'-k')
	# fill between the curves and add legend
	x=alti_bass_norm
	y1=px_bass_norm
	#plt.fill_between(alti_bass_norm,px_bass_norm,px_network_norm,where = px_bass_norm >= px_network_norm,facecolor='blue',interpolate=True,alpha = 0.3)
	plt.fill_between(x,y1,color='blue',interpolate=True,alpha = 0.1)
	x=alt_network_norm
	y1=px_network_norm
	plt.fill_between(x,y1,color='green',interpolate=True,alpha = 0.1)
	x=val_norm
	y1=min_norm
	plt.fill_between(x,y1,color='white',interpolate=True)	
	hh3 = 0.3 * rdr3
	hh1 = 0.3 * rdr1
	yr = (px_network_norm[int(float(hh3))] - min_norm[int(float(hh1))])/2 + min_norm[int(float(hh1))]
	plt.text(0.3,yr,r'Ir',style='oblique')
	hh3 = 0.5 * rdr3
	hh2 = 0.5 * rdr2
	yb = (px_bass_norm[int(float(hh2))] - px_network_norm[int(float(hh3))])/2 + px_network_norm[int(float(hh3))]
	plt.text(0.5,yb,r'Ib',style='oblique')
	plt.legend(('Hb','Hn','Hr'),'upper right',shadow=True)
	plt.xlabel('Normalized area (for Hb) or length (for Hn and Hr)')
	plt.ylabel('Height, normalized to the basin relief')
	plt.title(Title_string)
	Fig_name="Hypsometric_curves_R"+str(k+1)+".pdf"                           #generate figure file name
	savefig(graphs_path+"/"+Fig_name)                                         #Save figure	
	print " "
# End of doloop
# close f2w	
f2w.close()


# Calcul linear regression
print "________________________________________________"
print "Linear regression..."
# Init graph
plt.clf()
slope, intercept, r_value, p_value, std_err = sp.stats.linregress(lnA,R1w)
print 'slope SR =', slope
print 'r_value squared =', r_value**2
print 'p_value =', p_value
print 'standard deviation =', std_err

line = "slope SR =" + str(slope) + "  \n"
f0w.write(line)
line = "r_value squared =" + str(r_value**2) + "  \n"
f0w.write(line)
line = "p_value =" + str(p_value) + "  \n"
f0w.write(line)
line = "standard deviation =" + str(std_err) + "  \n"
f0w.write(line)
line = "  \n"
f0w.write(line)

# Plot regression
Title_string="SR_calculation"	
line = slope*lnA+intercept
plt.plot(lnA,line,'r-',lnA,R1w,'ro')
#pyplot.legend(('Hb','Hn','Hr'),'upper right',shadow=True)
txtlegend="Slope SR = "+str(slope)
plt.text(min(lnA),max(R1w)-0.1,r'Slope SR = '+str(slope),style='italic')
plt.text(min(lnA),max(R1w)-0.15,r'r2 = '+str(r_value**2),style='italic')
plt.xlabel('ln(A)')
plt.ylabel('R1w')
plt.title(Title_string)
Fig_name="SR.pdf"                                           #generate figure file name
savefig(graphs_path+"/"+Fig_name)                           #Save figure

# Give age of the perturbation of the drainage
# Calculate it with the Demoulin, GRL 2012 paper
t = 0.009 * slope**(-4.)
print "Age of the basin perturbation = ",t, "Ma"
line = "Age of the basin perturbation = " + str(t) + "Ma \n" 
f0w.write(line)
line = "  \n"
f0w.write(line)
f0w.close()

# END of code
