######!/usr/bin/env python
# -*- coding: UTF-8 -*-

############################
# Test to automate the extraction of watershed caracteristics from dem
############################

# Load python modules
#import scipy as sp : Not in Python arcgis !
import numpy as np
import os
from os import path, access, R_OK, mkdir         # W_OK for write permission
import sys
try:
	import Tkinter                                # To get info window
	import tkMessageBox, tkSimpleDialog           # To get info window
	from tkFileDialog import askopenfilename
except ImportError:
	# if Tkinter is not present, set textmode to 1 (no window popup)
	print "No Tkinter module. \n \n Please install it if you want to use graphic mode"
	textmode = 1
	pass


# Choose if Arcpy or Grass modules
# If Arcpy
import arcpy
from arcpy import env
from arcpy.sa import *
arcpy.CheckOutExtension("Spatial")
# check if arcgis 10.x


# If Grass




###############################################################################
###############################################################################
# PARMATERS TO BE CHANGED BY USER :
synthetic = "true"    # To define if the DEM is synthetic (= TRUE) or not (= FALSE)
factor = 20000      # For raster Calculator for f_Acc. HAS maybe to be calculated !
###############################################################################
###############################################################################
####### SET your ENVIRONNEMENT  #################
# set the working environnment
# You need to change the name of the variable to fit to your working environment
# Path where everything is stored
env.workspace = "//psf/Home/Documents/CODES-backup/ArcGIS/Demoulin-methodo/Essais/Xav-1R/SRSR_essai"
# Name of the input raster:
# The input raster must be projected in meters and in a projection that conserve the areas (e.g. Lambert, UTM,...) 		
input_raster = "fast"
# Name of the rasters used for the calculations
# (raster_fnme+"_corr", raster_fnme+"_fdir", raster_fnme+"_facc", raster_fnme+"_alti", raster_fnme+"_cost")
raster_fnme = "fast"
###############################################################################
###############################################################################

# Name of the generic river folders (will be increased by 1 for each next river) 
river_folder = env.workspace+"/"+"RProfils/R" 
raster_folder = env.workspace+"/"+"Rasters/" 
tmpf = env.workspace+"/TMP/"
# Load DEM
dem_fnme = raster_fnme
# Define folder where shp files will be stored
shp_fnme = env.workspace+"/"+"Shp_files/"+raster_fnme
# Check if Shp_files folder exists
if path.exists("Shp_files") == False:
	os.mkdir("Shp_files")
	# ERREUR DE DROITS D'aCCES
	
###################################################
# Function to find UTM zone from Lat Long
# Lat Long-UTM, UTM-Lat Long conversions
def LLtoUTM(Long, Lat, ReferenceEllipsoid = 23, zone = None):
	"""
	Converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
	East Longitudes are positive, West longitudes are negative
	North latitudes are positive, South latitudes are negative
	Lat and Long are in decimal degrees
	And This routine determines the correct UTM letter designator for the given latitude
	returns 'Z' if latitude is outside the UTM limits of 84N to 80S
	Written by Chuck Gantz- chuck.gantz@globalstar.com
	http://code.google.com/p/pyproj/issues/attachmentText?id=27&aid=-808841747718175643&name=UTM.py&token=cd8b549d525cfffddb3ade2d859ba68e
	"""
	
	from math import pi, sin, cos, tan, sqrt, radians, degrees
	
	_EquatorialRadius = 2
	_eccentricitySquared = 3
	
	_ellipsoid = [
	# id, Ellipsoid name, Equatorial Radius, square of eccentricity
	# first once is a placeholder only, To allow array indices to match id numbers
		[ 0, "Placeholder", 0.0, 0.0],
		[ 1, "European Datum 1950",    6378388.0000, 0.006723000],
		[ 2, "Airy",                   6377563.3960, 0.006670540],
		[ 3, "Australian National",    6378160.0000, 0.006694542],
		[ 4, "Bessel 1841",            6377397.1550, 0.006674372],
		[ 5, "Clarke 1866",            6378206.4000, 0.006768658],
		[ 6, "Clarke 1880",            6378249.1450, 0.006803511],
		[ 7, "Everest",                6377276.3452, 0.006637847],
		[ 8, "Fischer 1960 (Mercury)", 6378166.0000, 0.006693422],
		[ 9, "Fischer 1968",           6378150.0000, 0.006693422],
		[ 10, "GRS 1967",              6378160.0000, 0.006694605],
		[ 11, "GRS 1980",              6378137.0000, 0.006694380],
		[ 12, "Helmert 1906",          6378200.0000, 0.006693422],
		[ 13, "Hough",                 6378270.0000, 0.006722670],
		[ 14, "International",         6378388.0000, 0.006722670],
		[ 15, "Krassovsky",            6378245.0000, 0.006693422],
		[ 16, "Modified Airy",         6377340.1890, 0.006670540],
		[ 17, "Modified Everest",      6377304.0000, 0.006637847],
		[ 18, "Modified Fischer 1960", 6378155.0000, 0.006693422],
		[ 19, "South American 1969",   6378160.0000, 0.006694542],
		[ 20, "WGS 60",                6378165.0000, 0.006693422],
		[ 21, "WGS 66",                6378145.0000, 0.006694542],
		[ 22, "WGS 72",                6378135.0000, 0.006694318],
		[ 23, "WGS 84",                6378137.0000, 0.006694380]
	]
	
	# Reference ellipsoids derived from Peter H. Dana's website-
	# http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
	# Department of Geography, University of Texas at Austin
	# Internet: pdana@mail.utexas.edu
	# 3/22/95
	# Source
	# Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
	# 1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency
	
	a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
	eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
	k0 = 0.9996
	
	# Make sure the longitude is between -180.00 and 179.9
	LongTemp = (Long+180)-int((Long+180)/360)*360-180 #-180.00 .. 179.9
	
	LatRad = radians(Lat)
	LongRad = radians(LongTemp)
	
	if zone is None:
		ZoneNumber = int((LongTemp+180)/6)+1
	else:
		ZoneNumber = zone
	
	if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
		ZoneNumber = 32
	
	# Special zones for Svalbard
	if Lat >= 72.0 and Lat < 84.0:
		if  LongTemp >= 0.0  and LongTemp <  9.0:ZoneNumber = 31
		elif LongTemp >= 9.0  and LongTemp < 21.0: ZoneNumber = 33
		elif LongTemp >= 21.0 and LongTemp < 33.0: ZoneNumber = 35
		elif LongTemp >= 33.0 and LongTemp < 42.0: ZoneNumber = 37
		
	LongOrigin = (ZoneNumber-1)*6-180+3 #+3 puts origin in middle of zone
	LongOriginRad = radians(LongOrigin)
	
	# compute the UTM Zone from the latitude and longitude
	UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))
	
	eccPrimeSquared = (eccSquared)/(1-eccSquared)
	N = a/sqrt(1-eccSquared*sin(LatRad)*sin(LatRad))
	T = tan(LatRad)*tan(LatRad)
	C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
	A = cos(LatRad)*(LongRad-LongOriginRad)
	
	M = a*(
		(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256)*LatRad
		-(3*eccSquared/8+3*eccSquared*eccSquared/32+45*eccSquared*eccSquared*eccSquared/1024)*sin(2*LatRad)
		+(15*eccSquared*eccSquared/256+45*eccSquared*eccSquared*eccSquared/1024)*sin(4*LatRad)
		-(35*eccSquared*eccSquared*eccSquared/3072)*sin(6*LatRad)
	)
	
	UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
			+(5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
			+500000.0
	)
	
	UTMNorthing = (k0*(M+N*tan(LatRad)*(A*A/2
		+(5-T+9*C+4*C*C)*A*A*A*A/24
		+(61-58*T+T*T+600*C
		-330*eccPrimeSquared
		)*A*A*A*A*A*A/720)))
		
	if Lat < 0:
		UTMNorthing = UTMNorthing+10000000.0; #10000000 meter offset for southern hemisphere
    
	if   84 >= Lat >= 72: UTMLetter = 'X'
	elif 72 >  Lat >= 64: UTMLetter = 'W'
	elif 64 >  Lat >= 56: UTMLetter = 'V'
	elif 56 >  Lat >= 48: UTMLetter = 'U'
	elif 48 >  Lat >= 40: UTMLetter = 'T'
	elif 40 >  Lat >= 32: UTMLetter = 'S'
	elif 32 >  Lat >= 24: UTMLetter = 'R'
	elif 24 >  Lat >= 16: UTMLetter = 'Q'
	elif 16 >  Lat >= 8:  UTMLetter = 'P'
	elif  8 >  Lat >= 0:  UTMLetter = 'N'
	elif  0 >  Lat >=-8:  UTMLetter = 'M'
	elif -8 >  Lat >=-16: UTMLetter = 'L'
	elif -16 > Lat >=-24: UTMLetter = 'K'
	elif -24 > Lat >=-32: UTMLetter = 'J'
	elif -32 > Lat >=-40: UTMLetter = 'H'
	elif -40 > Lat >=-48: UTMLetter = 'G'
	elif -48 > Lat >=-56: UTMLetter = 'F'
	elif -56 > Lat >=-64: UTMLetter = 'E'
	elif -64 > Lat >=-72: UTMLetter = 'D'
	elif -72 > Lat >=-80: UTMLetter = 'C'
	else:                 UTMLetter = 'Z' # if the Latitude is outside the UTM limits
	
	return (UTMZone, UTMEasting, UTMNorthing, UTMLetter)
################################################


def check_raster(inraster,factor,operation):
	"""
	Function to check if the strahler order chosen is 
	"""
	is_valid=0
	print("Factor used to extract drainage is "+str(int(factor))+" \n Check if this is good")
	# Do a robust choice.
	while not is_valid :
		try:
			choice = int(aw_input("If the choice is good, enter 0 \n \
			                      If the choice is not good enough, \
			                      give a new value (Integer) for the factor"))
			is_valid = 1
		except ValueError, e :
			print("'%s' is not an integer!" % e.args[0].split(": ")[1])
	if choice != 0:
		factor = choice
		outsetnull = SetNull(inraster, inraster, "VALUE"+operation+str(factor))
		# Recusive function
		factor, outsetnull = check_raster(inraster, factor, operation)

	return factor,outsetnull
#################################################
#        END OF FUNCTION DEFINITIONS            #
#################################################



# Try if DEM is projected
# Get description
#############
print("###########################################################")
print("")
print ("Test raster projection...")
desc = arcpy.Describe(dem_fnme)
# Get spatial reference
sr = desc.spatialReference
# Else, project it in UTM format.
if sr.type != "Projected" and synthetic == false:
	print ("The DEM is not projected and not a synthetic dem. \n \
	        I project it in UTM, assuming that it is in Lat/Long - WSG84 reference frame")
	# Project it in UTM !! I Assume that of the raster is not projectec, that means that 
	# it is in geographical Lat-Long coordinates
	# read lat-long from the ll corner	
	lat1 = arcpy.GetRasterProperties_management(dem_fnme,"BOTTOM")
	lat2 = arcpy.GetRasterProperties_management(dem_fnme,"TOP")
	long1 = arcpy.GetRasterProperties_management(dem_fnme,"RIGHT")
	long2 = arcpy.GetRasterProperties_management(dem_fnme,"LEFT")
	# Calcul center of the raster
	lat = lat1 + abs(lat2 - Lat1)/2
	long = long1 + abs(long2 - Long1)/2
	# find Hemisphere
	if lat < 0: 
		hemisphere = 'S'	
	else:
		hemisphere = 'N'		
	UTMZone, UTMEasting, UTMNorthing, UTMLetter = LLtoUTM(long, lat)
	outCS = arcpy.SpatialReference('NAD 1983 UTM Zone '+UTMzone+hemisphere)	
	try:
		arcpy.ProjectRaster_managment(dem_fnme, dem_fnme_proj, outCS)	      
		dem_fnme = dem_fnme_proj
		# In that case, we should build (or rebuild) all the existing rasters.
		# Thus, clean the whole workspace
		os.rmdir(raster_folder)
		os.mkdir(raster_folder)
		os.rmdir(env.workspace+"/RProfils")
		os.mkdir(env.workspace+"/RProfils")
		os.rmdir(tmpf)
		os.mkdir(tmpf)
	except:
		print "Project Raster failed"
		print arcpy.GetMessages()
print("Main raster is projected")
print("")
print("###########################################################")
# Check if working folder is clean or not. If this is not clean, it is not possible to rewrite on it, and thus the automatic code cannot be ran
print("Check if working folder is clean or not. \n If this is not clean, it is not possible to rewrite on it, \n and thus the automatic code cannot be ran")
if path.isfile(raster_folder+"in_stream_r") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(raster_folder+"in_stream_r")))    # if yes, print error message and exit
if path.isfile(raster_folder+"stream_order") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(raster_folder+"stream_order")))    # if yes, print error message and exit
if path.isfile(shp_fnme+'_streams.shp') == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+'_streams.shp')))    # if yes, print error message and exit
if path.isfile(shp_fnme+"_strorder.shp") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+"_strorder.shp")))    # if yes, print error message and exit
if path.isfile(shp_fnme+"_strorder2.shp") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+"_strorder2.shp")))    # if yes, print error message and exit
if path.isfile(shp_fnme+"_intersect2.shp") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+"_intersect2.shp")))    # if yes, print error message and exit
if path.isfile(shp_fnme+"_intersect.shp") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+"_intersect.shp")))    # if yes, print error message and exit
if path.isfile(shp_fnme+"_outlets.shp") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(shp_fnme+"_outlets.shp")))    # if yes, print error message and exit
if path.isfile(env.workspace+"/River-caract+txt") == True:
	sys.exit('File {FileNa} exists \n Clean your workspace'.format(FileNa=str(env.workspace+"/River-caract+txt")))    # if yes, print error message and exit

# Check if basics rasters exists
print("Check if flow direction, accumulation, alti and cost rasters exist... \n \
       If there is a need to build them, be patient...")
if path.isfile(raster_folder+raster_fnme+"_corr") == False or access(raster_folder+raster_fnme+"_corr", R_OK) == False:
	print("It seems that rasters files are not prepared ! \n We are doing it...")
	if path.isfile(raster_folder+raster_fnme+"_corr") == True and access(raster_folder+raster_fnme+"_corr", R_OK) == False:
		os.chmod(raster_folder+raster_fnme+"_corr",711)
	else:
		#sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_corr")))    # if yes, print error message and exit
		print "      extract values > 0 m"
		# Con
		outcon = Con(raster_folder+input_raster,raster_folder+input_raster,"VALUE" >= 0)
		print("      Fill holes")
		# Fill holes
		outfill = Fill(outcon)
		outfill.save(raster_folder+raster_fnme+"_corr")
if path.isfile(raster_folder+raster_fnme+"_fdir") == False or access(raster_folder+raster_fnme+"_fdir", R_OK) == False:
	#sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_fdir")))    # if yes, print error message and exit
	if path.isfile(raster_folder+raster_fnme+"_fdir") == True and access(raster_folder+raster_fnme+"_fdir", R_OK) == False:
		os.chmod(raster_folder+raster_fnme+"_fdir",711)
	else:
		print("      Flow direction")
		# Flow direction
		outflowdir = FlowDirection(raster_folder+raster_fnme+"_corr")
		outflowdir.save(raster_folder+raster_fnme+"_fdir")
if path.isfile(raster_folder+raster_fnme+"_facc") == False or access(raster_folder+raster_fnme+"_facc", R_OK) == False:
	#sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_facc")))    # if yes, print error message and exit
	if path.isfile(raster_folder+raster_fnme+"_facc") == True and access(raster_folder+raster_fnme+"_facc", R_OK) == False:
		os.chmod(raster_fnme+"_facc",711)
	else:
		print("      Flow accumulation")
		# Flow accumulation
		outflowacc = FlowAccumulation(raster_folder+raster_fnme+"_fdir","","INTEGER")
		outflowacc.save(raster_folder+raster_fnme+"_facc")
if path.isfile(raster_folder+raster_fnme+"_alti") == False or access(raster_folder+raster_fnme+"_alti", R_OK) == False:
	#sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_alti")))    # if yes, print error message and exit
	if path.isfile(raster_folder+raster_fnme+"_alti") == True and access(raster_folder+raster_fnme+"_alti", R_OK) == False:
		os.chmod(raster_folder+raster_fnme+"_alti",711)
	else:
		print("      Alti Raster")
		# Alti raster
		# 62 correspond to the minimum size of a basin at 0.5 km2.
		outcon = Con(raster_folder+raster_fnme+"_facc",raster_folder+raster_fnme+"_corr","VALUE" >= 62)
		outcon.save(raster_folder+raster_fnme+"_alti")
if path.isfile(raster_folder+raster_fnme+"_cost") == False or access(raster_folder+raster_fnme+"_cost", R_OK) == False:
	#sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_cost")))    # if yes, print error message and exit
	if path.isfile(raster_folder+raster_fnme+"_cost") == True and access(raster_folder+raster_fnme+"_cost", R_OK) == False:
		os.chmod(raster_folder+raster_fnme+"_cost",711)
	else:
		print("      Cost Raster")
		# Cost raster
		outcon = Con(raster_folder+raster_fnme+"_facc",1,"VALUE" >= 62)
		outcon.save(raster_folder+raster_fnme+"_cost")
		print(" ")
print("flow direction, accumulation, alti and cost rasters exist !")
print(" ")

#### See http://gis4geomorphology.com/watershed/ for the description of the different steps ####
print("###########################################################")
print("Begin streams extraction procedure...")
print(" ")
print("Compute the stream network...")
# Calcul the threshold to extract drainage
factor = arcpy.GetRasterProperties_management(raster_folder+raster_fnme+"_facc", "MAXIMUM").getOutput(0)
factor = factor/50
outsetnull = SetNull(raster_folder+raster_fnme+"_facc",raster_folder+raster_fnme+"_facc","VALUE > "+str(factor))
# check if this is OK
factor, outsetnull = check_raster(raster_folder+raster_fnme+"_facc", factor, " > ")
outsetnull.save(raster_folder+"in_stream_r")
# Convert Stream Pixels to Shapefile
StreamToFeature (raster_folder+"in_stream_r", raster_folder+raster_fnme+'_fdir', shp_fnme+'_streams.shp', "NO_SIMPLIFY")	


# Delimitate several watersheds:
# Run StreamLink tool
##arcpy.StreamLink_sa(path+"in_stream_raster", path+"_fdir", Output_raster_3)
#outstreamlink = arcpy.StreamLink("in_stream_raster2", raster_fnme+"_fdir")
#outstreamlink.save(raster_fnme+"_streamlink")

print("")
print("Compute Strahler order")
# Find Stream orders
orderMethod = "STRAHLER"
outstreamorder = StreamOrder(raster_folder+"in_stream_raster", raster_folder+raster_fnme+'_fdir', orderMethod)
outstreamorder.save(raster_folder+"stream_order")
# convert to shapefile ????

print("")
print("Extract outlets of interest")
# Find outlet of each relevent stream 
#1) Convert stream order raster to shp for stream order = strahler_threshold
strahler_threshold = 3
outsetnull = SetNull(raster_folder+"stream_order",raster_folder+"stream_order","VALUE <> "+str(strahler_threshold))
# Check to know if the strahler order choosen is good or no !!!
strahler_threshold, outsetnull = check_raster(raster_folder+"stream_order",strahler_threshold," <> ")
StreamToFeature(outsetnull, outsetnull,shp_fnme+"_strorder.shp", "SIMPLIFY")

# Find rivers with stream_order = strahler_threshold+1
outsetnull = SetNull(raster_folder+"stream_order",raster_folder+"stream_order","VALUE <> "+str(strahler_threshold+1))
StreamToFeature(outsetnull, outsetnull,shp_fnme+"_strorder2.shp", "SIMPLIFY")
# Find intersections between stream_order_n rivers and the boundary of the working area = outlets
arcpy.Intersect_analysis([shp_fnme+"_strorder2.shp",shp_fnme+"_strorder.shp"],shp_fnme+"_intersect2.shp","ALL",40,"POINT")
# Build contour of the working area
outcon = Con(raster_folder+raster_fnme+"_corr",1,raster_folder+raster_fnme+"_corr" > 0)
arcpy.RasterToPolygon_conversion(outcon,shp_fnme+"_rextent")
# Find intersections between stream_order_n rivers and the boundary of the working area = outlets
arcpy.Intersect_analysis([shp_fnme+"_rextent.shp",shp_fnme+"_strorder.shp"],shp_fnme+"_intersect.shp","ALL",40,"POINT")
# Merge the two intersect shapefiles
arcpy.Merge_management([shp_fnme+"_intersect.shp",shp_fnme+"_intersect2.shp"],shp_fnme+"_outlets.shp")
# remove intersects_files
arcpy.Delete_management(shp_fnme+"_intersect.shp")
arcpy.Delete_management(shp_fnme+"_intersect2.shp")
# Clean the shapefile by removing double points
arcpy.AddXY_management(shp_fnme+"_outlets.shp")
arcpy.DeleteIdentical_management(shp_fnme+"_outlets.shp",["POINT_X","POINT_Y"])

river = open(env.workspace+"/River-caract+txt", "w")
river.write("river_nb \t L_watershed \n")
#3) split Stream_X into shp of 1 point each, and store it into a folder Ri
nbriver = int(arcpy.GetCount_management(shp_fnme+"_outlets.shp").getOutput(0))
print("Working on "+ str(nbriver)+" basins")
for i in range (0,nbriver):
	if path.exists(tmpf) != False:
		os.rmdir(tmpf)
	os.mkdir(tmpf)	

	print("         Outlet of basin "+str(i+1))
	path_river = env.workspace+"/RProfils/R"+str(i+1)+"/shp_files/"
	path_river_r = env.workspace+"/RProfils/R"+str(i+1)+"/rasters/"
	os.mkdir(env.workspace+"/RProfils/R"+str(i+1))
	os.mkdir(path_river)
	os.mkdir(path_river_r)
	os.mkdir(path_river+"/tmp")
	arcpy.FeatureClassToFeatureClass_conversion(shp_fnme+"_outlets.shp", path_river, \
	      "exu_R"+str(i+1)+".shp","FID = "+str(i))
	# Extract watershed/basins of each relevant stream
	print("         build basin "+str(i+1))
	outsnappour_exu = SnapPourPoint(path_river+"exu_R"+str(i+1)+".shp", \
	                  raster_folder+raster_fnme+"_facc",100,"Id")
	outsnappour_exu.save(path_river+"/tmp/g_exu_R"+str(i+1))
	outwatershed = Watershed(raster_folder+raster_fnme+"_fdir", outsnappour_exu)
	outwatershed.save(path_river_r+raster_fnme+"_wshed"+str(i+1))
	arcpy.RasterToPolygon_conversion(path_river_r+raster_fnme+"_wshed"+str(i+1), \
	                                 path_river_r+raster_fnme+"_wshed"+str(i+1)+".shp", "NO_SIMPLIFY")	
	                              	
	print("         Extract source point "+str(i+1))
	# Extract Source/head channel points of relevant watersheds
	# mask the flow dir raster with the basin geometry
	outextractmask = ExtractByMask(raster_folder+raster_fnme+"_facc", path_river_r+raster_fnme+"_wshed"+str(i+1)+".shp") 
	## Find the cell with the highest flow length
	#outflowlength = FlowLength(outextractmask, "DOWNSTREAM") # Calculer le In_flow_Raster pour chaque bassin
	## get X/Y/flowlength of the cell with the max
	#maxflowLobject = arcpy.GetRasterProperties_management(outflowlength, "MAXIMUM")
	#maxflow = maxflowLobject.getOutput(0)
	#outsetnull = SetNull(outflowlength, 1, "VALUE < "+str(maxflow))
	#arcpy.RasterToPoint_conversion(outsetnull, path_river+"src_R"+str(i+1)+".shp")
	## PROBLEM? IT DOES NOT WORK....

	# BE CAREFUll : the value threshold is set to 250, but it can maybe less.
	#               The idea is to find the closest point of the basin boundary !
	outsetnull = SetNull(outextractmask,outextractmask,"VALUE < 250")
	StreamToFeature(outsetnull, outsetnull,tmpf+"tmp.shp", "SIMPLIFY")
	arcpy.FeatureVerticesToPoints_management(tmpf+"tmp.shp", tmpf+"tmp2.shp","END")
	arcpy.PointDistance_analysis(tmpf+"tmp2.shp", path_river+"exu_R"+str(i+1)+".shp", \
		                             env.workspace+"RProfils/R"+str(i+1)+"table_acc_dist"+str(i+1))
	# Find FID of the MAX_DISTANCE,
	arcpy.Statistics_analysis(env.workspace+"/RProfils/R"+str(i+1)+"/table_acc_dist"+str(i+1)+".dbf", \
	                          tmpf+"stats1", [["DISTANCE", "MAX"]])
	rows = arcpy.SearchCursor(tmpf+"stats1")
	row = rows.next()
	maxd = row.MAX_DISTANCE
	arcpy.TableSelect_analysis(tmpf+"stats1", tmpf+"select1", "DISTANCE = "+str(maxd))
	rows = arcpy.SearchCursor(tmpf+"select1")
	row = rows.next()
	fidr = row.INPUT_FID
	# Extract point with the FID of the MAX_DISTANCE !
	arcpy.FeatureClassToFeatureClass_conversion(tmpf+"tmp2.shp", path_river, "src_R"+str(i+1)+".shp","FID = "+str(fidr))
	
	print("         Calcul length of the basin "+str(i+1))      
	# Find max length of the basin (not the length of the stream): for each watershed, 
	#                compute the distance between the outlet and each cell of the 
	#                watershed limit. Take the longest
	arcpy.PolygonToLine_management(path_river_r+raster_fnme+"_wshed"+str(i+1)+".shp", \
	            tmpf+raster_fnme+"+_wshed_line"+str(i+1)+".shp")
	arcpy.FeatureVerticesToPoints_management(tmpf+raster_fnme+"_wshed_line"+str(i+1)+".shp", \
	            tmpf+raster_fnme+"_wshed_pts"+str(i+1)+".shp","END")
	            # It seems to be correct, but it gives just one point,
	            # the line is not cut in several lines
	            # Maybe, we have to check that this is always correct !
	arcpy.PointDistance_analysis(tmpf+raster_fnme+"_wshed_pts"+str(i+1)+".shp", \
		                             path_river+"exu_R"+str(i+1)+".shp", \
		                             tmpf+"table_dwshed"+str(i+1))
	arcpy.Statistics_analysis(tmpf+"table_dwshed"+str(i+1)+".dbf", \
	                          tmpf+"stats2", [["DISTANCE", "MAX"]])
	rows = arcpy.SearchCursor(tmpf+"stats2")
	row = rows.next()
	Lwatershedmax = row.MAX_DISTANCE
	# clean
	print("")
	print("Cleaning.............")
	os.rmdir(tmpf)
	print("         write in file the line "+str(i+1))
	# write the line in file River-caract.txt corresponding to the watershed i+1
	river.write(str(i+1)+"\t"+str(Lwatershedmax)+"\n")	
river.close()

print("****************************")
print("exu_R and src_R Shape files have been written for each river of interest")
print("File River-caract.txt writen \n"+str(nbriver)+" rivers have been selected")
print("****************************")

	

