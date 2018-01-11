######!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

###############################################################################
#                                                                             #
# ArcGIS script to extract data to apply Demoulin approach on rivers profiles #
#      Based on the work of Demoulin, 2011, 2012, and Thomas Croissant        #
#                              Work on ArcGIS 10.1                            #
#                 modified by X. Robert, Grenoble, 08/04/2013                 #
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
#		>>> execfile("2-script-ArcGIS.py")


#import system modules:
import arcpyfrom arcpy import env
from arcpy.sa import *
import os
from os import path, access, R_OK, mkdir         # W_OK for write permission
import sys

# TO DO : 
#	- If some rivers have been already done, it crashs. So 2 solutions :
#		* We do a check if it has already been done, and if yes,
#		  we erase the folder tmp and the dbf files
#		* We do a check if it has already been done, and if yes,
#         we jump to the next river (quicker, but do not work if it has previously
#         crashed
#       * It is maybe possible to do a mix
#         (ex, check if a file is missing, and if yes, do not skeep, but erase)
 
###############################################################################
###############################################################################
####### SET your ENVIRONNEMENT  #################
# set the working environnment
# You need to change the name of the variable to fit to your working environment
		# Path where everything is stored
#env.workspace = "//psf/Home/Documents/SIG/Peru/Dem/Desmoulins-essai/SRTM"
env.workspace = "//psf/Home/Documents/CODES-backup/ArcGIS/Demoulin-methodo/Essais/Xav-1R/SRSR"		
		# Name of the rasters used (raster_fnme+"_corr", raster_fnme+"_fdir", raster_fnme+"_facc", raster_fnme+"_alti", raster_fnme+"_cost")
raster_fnme = "fast"
		# Name of the generic river folders (will be increased by 1 for each next river) 
river_folder = "RProfils/R"  
		# Name of the text-file with the rivers caracteristics
		#     first line =  header [river_nb,L_watershed]
		#     then, one line per river, with 2 columns [river_nb,L_watershed]
		#     tab separation between the columns
Rivercaract = "River-caract.txt"
		
######## END of modifications ###################
###############################################################################
###############################################################################
print "ArcGIS script to extract data to apply Demoulin approach"
print "on rivers profiles"
print "X. Robert, Grenoble, 06/2013"
print "  "
print "  "

print "Checking if paths and input files exists..."
# Folder where rasters are stored
raster_fnme ="Rasters/"+raster_fnme
# Check if workspace and input rasters exists
if path.exists(env.workspace) == False:            #Check if the path exist
	sys.exit('Path workspace {PathNa} does not exist'.format(PathNa=str(env.worspace)))    # if not, print error message and exit
# Check if raster files needed exist
#_corr _fdir _facc _alti _cost
if path.isfile(raster_fnme+"_corr") == False and access(raster_fnme+"_corr", R_OK) == False:
	sys.exit('File {FileNa} does not exist'.format(FileNa=str(raster_fnme+"_corr")))    # if not, print error message and exit
if path.isfile(raster_fnme+"_fdir") == False and access(raster_fnme+"_fdir", R_OK) == False:
	sys.exit('File {FileNa} does not exist'.format(FileNa=str(raster_fnme+"_fdir")))    # if not, print error message and exit
if path.isfile(raster_fnme+"_facc") == False and access(raster_fnme+"_facc", R_OK) == False:
	sys.exit('File {FileNa} does not exist'.format(FileNa=str(raster_fnme+"_facc")))    # if not, print error message and exit
if path.isfile(raster_fnme+"_alti") == False and access(raster_fnme+"_alti", R_OK) == False:
	sys.exit('File {FileNa} does not exist'.format(FileNa=str(raster_fnme+"_alti")))    # if not, print error message and exit
if path.isfile(raster_fnme+"_cost") == False and access(raster_fnme+"_cost", R_OK) == False:
	sys.exit('File {FileNa} does not exist'.format(FileNa=str(raster_fnme+"_cost")))    # if not, print error message and exit


# read the file River-charact.txt
# Should contain Rivers numbers with length of the watershed. See example provided
print ('Reading {FileNa}'.format(FileNa=str(Rivercaract)))
if path.isfile(env.workspace+"/"+Rivercaract) and access(env.workspace+"/"+Rivercaract, R_OK):            #Check if the file exist
	river = open(env.workspace+"/"+Rivercaract, "r")
	nd = 0                                                            #nd = number of rivers to read
	Lwatershed = []
	riversnb = []
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
	print "If you have lots of river, check the first iteration, and then grab a beer, it can be long !"
	print "  "
else:                                                                  # If file does not exist, print error message and quit
	print('ERROR : File {FileNa} does not exist'.format(FileNa=str(Rivercaract)))
	print('in the folder',env.workspace) 
	sys.exit('File {FileNa} should contain an hearder (first line), and then the number of the river and the length of the watershed'.format(FileNa=str(Rivercaract))) 

##################
# River Stepping #
##################
#nd = 0
for i in range(nd):
	print "Working on river", i+1
	River = env.workspace+"/"+river_folder+str(i+1)
	# Folder where shapefiles are stored
	shp_folder = River+"/shp_files"
	#Check if Folders exists
	if path.exists(River) == False:
		# If it does not exist (this one should exist), create it
		mkdir(River)
		print "Creating folder ",River
	if path.exists(River+"/tmp") == False:
		# If it does not exist, create it
		mkdir(River+"/tmp")
		print "Creating folder ",River+"/tmp"
			
	# Check if shapefiles needed exist
	# src_R exu_R
	if path.isfile(shp_folder+"/exu_R"+str(i+1)+".shp") == False and access(shp_folder+"/exu_R"+str(i+1)+".shp", R_OK) == False:
		sys.exit('File {FileNa} does not exist'.format(FileNa="exu_R"+str(i+1)+".shp"))    # if not, print error message and exit
	if path.isfile(shp_folder+"/src_R"+str(i+1)+".shp") == False and access(shp_folder+"/src_R.shp", R_OK) == False:
		sys.exit('File {FileNa} does not exist'.format(FileNa="src_R"+str(i+1)+".shp"))    # if not, print error message and exit

	print "   Snap Pour Point" 
	# SnapPourPoint : Snaps pour points to the cell of highest flow accumulation within a specified distance.
	outSnapPour = SnapPourPoint(shp_folder+"/exu_R"+str(i+1)+".shp",env.workspace+"/"+raster_fnme+"_facc",100,"Id")
	outSnapPour.save(River+"/tmp/g_exu_R"+str(i+1))

	# SnapPourPoint : Snaps pour points to the cell of highest flow accumulation within a specified distance.	outSnapPour = SnapPourPoint(shp_folder+"/src_R"+str(i+1)+".shp",env.workspace+"/"+raster_fnme+"_facc",100,"Id")
	outSnapPour.save(River+"/tmp/g_src_R"+str(i+1))

	# Set the extent environment using a keyword.
	arcpy.env.extent = "MAXOF"

	print "   Cost calculations"
	# Calculates the least accumulative cost distance for each cell to the nearest source over a cost surface	outCostDist = CostDistance(River+"/tmp/g_src_R"+str(i+1),env.workspace+"/"+raster_fnme+"_cost",200000)	outCostDist.save(River+"/tmp/R"+str(i+1)+"_costdi")
	# Set the extent environment using a keyword.
	arcpy.env.extent = "MAXOF"

	## Defines the neighbor that is the next cell on the least accumulative cost path to the nearest source.	outCostBackLink=CostBackLink(River+"/tmp/g_src_R"+str(i+1),env.workspace+"/"+raster_fnme+"_cost",100000)
	outCostBackLink.save(River+"/tmp/R"+str(i+1)+"_costbl")	# Calculates the least-cost path from a source to a destination.	outCostPath=CostPath(River+"/tmp/g_exu_R"+str(i+1),River+"/tmp/R"+str(i+1)+"_costdi",River+"/tmp/R"+str(i+1)+"_costbl")
	outCostPath.save(River+"/tmp/R"+str(i+1))

	print "   Extraction by masks"		# Extracts the cells of a raster that correspond to the areas defined by a mask.	outExtractByMask=ExtractByMask(env.workspace+"/"+raster_fnme+"_alti", River+"/tmp/R"+str(i+1))
	outExtractByMask.save(River+"/tmp/R"+str(i+1)+"_alti")
	outExtractByMask=ExtractByMask(River+"/tmp/R"+str(i+1)+"_costdi",River+"/tmp/R"+str(i+1))
	outExtractByMask.save(River+"/tmp/R"+str(i+1)+"_dist")
	outSingleOutputMapAlgebra=Int(Raster(River+"/tmp/R"+str(i+1)+"_dist")/200) + 1
	outSingleOutputMapAlgebra.save(River+"/tmp/R"+str(i+1)+"_segm")

	# Summarizes the values of a raster within the zones of another dataset and reports the results to a table.	outZonalStatisticsAsTable=ZonalStatisticsAsTable(River+"/tmp/R"+str(i+1)+"_segm","VALUE",River+"/tmp/R"+str(i+1)+"_alti",River+"/tmp/R"+str(i+1)+"_profil")

	print "   Watershed"
	# Determines the contributing area above a set of cells in a raster.	outWatershed=Watershed(raster_fnme+"_fdir", River+"/tmp/g_exu_R"+str(i+1))
	outWatershed.save(River+"/tmp/R"+str(i+1)+"_bass")
	outExtractByMask=ExtractByMask(env.workspace+"/"+raster_fnme+"_corr", River+"/tmp/R"+str(i+1)+"_bass")
	outExtractByMask.save(River+"/tmp/R"+str(i+1)+"_demb")
	outExtractByMask=ExtractByMask(env.workspace+"/"+raster_fnme+"_alti", River+"/tmp/R"+str(i+1)+"_bass")
	outExtractByMask.save(River+"/tmp/R"+str(i+1)+"_demr")

	print "   Write .dbf tables"
	# write the .dbf tables	arcpy.TableToTable_conversion(River+"/tmp/R"+str(i+1)+"_profil", River,"R"+str(i+1)+".dbf")
	# "VALUE VALUE true false false 4 Long 0 0 ,First,
	#,PATH\R1_profil,VALUE,-1,-1;COUNT COUNT true false false 4 Long 0 0 ,First,#,PATH\R1_profil,COUNT,-1,-1;MIN MIN true false false 4 Long 0 0 ,First,
	#,PATH\R1_profil,MIN,-1,-1;MAX MAX true false false 4 Long 0 0 ,First,
	#,PATH\R1_profil,MAX,-1,-1" 
	# PATH\profils\R1.dbf
	arcpy.TableToTable_conversion(River+"/tmp/R"+str(i+1)+"_demb", River, "R"+str(i+1)+"demb.dbf")
	# "VALUE VALUE true false false 4 Long 0 0 ,First,
	#,PATH\R1_demb,VALUE,-1,-1;COUNT COUNT true false false 4 Long 0 0 ,First,
	#,PATH\R1_demb,COUNT,-1,-1" 
	# PATH\profils\Rdemb.dbf
	arcpy.TableToTable_conversion(River+"/tmp/R"+str(i+1)+"_demr", River, "R"+str(i+1)+"demr.dbf")
	# "VALUE VALUE true false false 4 Long 0 0 ,First,
	#,PATH\R1_demr,VALUE,-1,-1;COUNT COUNT true false false 4 Long 0 0 ,First,
	#,PATH\R1_demr,COUNT,-1,-1" 
	# PATH\profils\Rdemr.dbf
	print " "

print "Script finished"# End of script
