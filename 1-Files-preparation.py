######!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

###############################################################################
#                                                                             #
#         ArcGIS script to prepare files to apply Demoulin approach           #
#                            on rivers profiles                               # 
#                            Work on ArcGIS 10.x                              #
#                       X. Robert, Grenoble, 06/2013                          #
#                                                                             #
###############################################################################

# You need to run this python script in ArcGIS (python terminal).
# You need to have enabled the extension "spatial analyst"
#    This script does NOT check it
#
# Follow the structure of the data base needed
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
# you need to go in the working folder
#        >>> import os
#        >>> os.chdir("your/path")
#
# And then, to run the script
#		>>> execfile("1-Files-preparation.py")

#import system modules:
import arcpy                                     # needs ARCGIS on windows
from arcpy import env
from arcpy.sa import *
import os
from os import path, access, R_OK, mkdir         # W_OK for write permission
import sys

###############################################################################
###############################################################################
####### SET your ENVIRONNEMENT  #################
# set the working environnment
# You need to change the name of the variable to fit to your working environment
# Path where everything is stored
env.workspace = "//psf/Home/Documents/CODES-backup/ArcGIS/Demoulin-methodo/Essais/Xav-1R/SRSR"
# Name of the input raster:
# The input raster must be projected in meters and in a projection that conserve the areas (e.g. Lambert, UTM,...) 		
input_raster = "fast"
# Name of the rasters used for the calculations
# (raster_fnme+"_corr", raster_fnme+"_fdir", raster_fnme+"_facc", raster_fnme+"_alti", raster_fnme+"_cost")
raster_fnme = "fast"
# Name of the generic river folders (will be increased by 1 for each next river) 
river_folder = "RProfils/R"  
		
######## END of modifications ###################
###############################################################################
###############################################################################
# Folder where will be stored the graphs (better not to change it)
graphs_path = "Graphs" 
	
print "ArcGIS script to prepare files to apply Demoulin approach"
print "on rivers profiles"
print "X. Robert, Grenoble, 06/2013"
print "  "
print "  "

# Check and build the structure if not exists
# Check if workspace and input rasters exists
print "Check existing struture..."
if path.exists(env.workspace) == False:            #Check if the path exist
	sys.exit('Path workspace {PathNa} does not exist'.format(PathNa=str(env.workspace)))    # if not, print error message and exit
# Check if the structure of the workspace exists. If not, create it.
if path.exists(env.workspace+"/Rasters") == False:
	os.mkdir(env.workspace+"/Rasters")
if path.exists(env.workspace+"/"+river_folder[-2]) == False:
	os.mkdir(env.workspace+"/"+river_folder[-2])
#### WHy the next 6 lines ???? --> If not automatic, this is to give the folders structure !
#if path.exists(env.workspace+"/"+river_folder+"1") == False:
#	os.mkdir(env.workspace+"/"+river_folder+"1")
#if path.exists(env.workspace+"/"+river_folder+"1/shp_files") == False:
#	os.mkdir(env.workspace+"/"+river_folder+"1/shp_files")
#if path.exists(env.workspace+"/"+river_folder+"1/rasters") == False:
#	os.mkdir(env.workspace+"/"+river_folder+"1/rasters")
	
if path.exists(env.workspace+"/Shp_files") == False:
	os.mkdir(env.workspace+"/Shp_files")
if path.exists(env.workspace+"/ks_auto") == False:
	os.mkdir(env.workspace+"/ks_auto")
if path.exists(env.workspace+"/ks_auto/Raw_outputs") == False:
	os.mkdir(env.workspace+"/ks_auto/Raw_outputs")
if path.exists(env.workspace+"/ks_auto/shp_files") == False:
	os.mkdir(env.workspace+"/ks_auto/shp/files")
# Verify if the folder for Graphs exists
if path.exists(env.workspace+"/"+graphs_path) == False:
	# If it do not exists, create it
	mkdir(env.workspace+"/"+graphs_path)
if path.exists(env.workspace+"/TMP") == False:	
	os.mkdir("TMP")
	
raster_fnme = env.workspace+"/Rasters/"+raster_fnme

# Check if raster files exist
#_corr _fdir _facc _alti _cost
if path.isfile(raster_fnme+"_corr") == True and access(raster_fnme+"_corr", R_OK) == True:
	sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_corr")))    # if yes, print error message and exit
if path.isfile(raster_fnme+"_fdir") == True and access(raster_fnme+"_fdir", R_OK) == True:
	sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_fdir")))    # if yes, print error message and exit
if path.isfile(raster_fnme+"_facc") == True and access(raster_fnme+"_facc", R_OK) == True:
	sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_facc")))    # if yes, print error message and exit
if path.isfile(raster_fnme+"_alti") == True and access(raster_fnme+"_alti", R_OK) == True:
	sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_alti")))    # if yes, print error message and exit
if path.isfile(raster_fnme+"_cost") == True and access(raster_fnme+"_cost", R_OK) == True:
	sys.exit('File {FileNa} exists'.format(FileNa=str(raster_fnme+"_cost")))    # if yes, print error message and exit

# Check if the raster is projected/synthetic or not. If not, we have to project it !!!
# See file Essai Automatic for a try !




# Prepare the different files
print "extract values > 0 m"
# Con
outcon = Con(input_raster,input_raster,"VALUE" >= 0)
outcon.save(raster_fnme)
################# DO not work
# make clean
#arcpy.CopyRaster_management(input_raster,"Rasters/"+input_raster+"_input")
#print 'Original Raster {FileNa} have been copied to Rasters/.'.format(FileNa=str(raster_fnme))
#print 'You can delete the Raster {FileNa} from the parent folder'.format(FileNa=str(raster_fnme))
################# DO not work


print(" ")
print("Fill holes")
# Fill holes
outfill = Fill(raster_fnme)
outfill.save(raster_fnme+"_corr")
print("Flow direction")
# Flow direction
outflowdir = FlowDirection(raster_fnme+"_corr")
outflowdir.save(raster_fnme+"_fdir")
print("Flow accumulation")
# Flow accumulation
outflowacc = FlowAccumulation(raster_fnme+"_fdir","","INTEGER")
outflowacc.save(raster_fnme+"_facc")
print("Alti Raster")
# Alti raster
# 62 correspond to the minimum size of a basin at 0.5 km2.
outcon = Con(raster_fnme+"_facc",raster_fnme+"_corr","VALUE" >= 62)
outcon.save(raster_fnme+"_alti")
print("Cost Raster")
# Cost raster
outcon = Con(raster_fnme+"_facc",1,"VALUE" >= 62)
outcon.save(raster_fnme+"_cost")
print(" ")
print("Create the shapefiles in the first river" )
# Create the shapefiles in the first river
arcpy.CreateFeatureclass_management("RProfils/R1/shp_files", "src_R1.shp", "POINT")
arcpy.CreateFeatureclass_management("RProfils/R1/shp_files", "exu_R1.shp", "POINT")

print(" ")
print("Script finished")
print(" ")
print("You now have to choose your rivers,")
print("mesure their watershed length (create 'River_caract.txt'),")
print("and to build the 'src_R.shp' and 'exu_R.shp'")
print("See 'README.txt' file")

# END of Script