ArcGIS script to extract data to apply Demoulin approach on rivers profiles
===========================================================================

Based on the work of Demoulin, 2011, 2012, and Thomas Croissant        
                              
It works on ArcGIS 10.x                            

X. Robert, Grenoble, 06/2013                        


This tool have been set up to automatize the calcul of the Demoulin indices from a DEM raster. It has 3 main files:

1. files-preparation.py: create the structure of the workspace (if not already created)
						create the five rasters needed for the next scripts							

2. Scrip-ArcGIS.py: Create the .dbf tables needed for the calculs of indices

3. Calculindices_Demoulin.py: Calculate the indices
	
They have to be ran in the order given by the numbers


REQUIREMENTS:
-------------
	
	- ArcGIS 10.x with python 2.7 and spatial analyst toolbox
	
	- Additional Python 2.7 independent from the one provided with ArcGIS (I didn't manage to install the required python packages on the ArcGIS Python)
	
	- package numpy 1.7 or higher (install or update via easy_install)
	
	- package Scipy 0.12 or higher
	
	- package dbfpy (need to be installed from http://dbfpy.sourceforge.net, see python package installation tutorial)
	
	- (install or update via easy_install)


How to run the scripts: Cookbook
--------------------------------

1. Before to run the script 1-files-preparation.py, you need a DEM raster projected in a projection in meters, with area conservations (e.g. UTM). To project a raster, use the ArcGIS tool Data_Managment_Tools --> Raster --> Project_Raster

2. In ArcGIS, Open the Python terminal:	you need to go in the working folder

	.. code-block:: python
    
		>>> import os

		>>> os.chdir("your/path")

	Edit the script 1-Files-preparation.py (with a text editor)
	
	- Go to ####### SET your ENVIRONNEMENT  #################
	
	- Give the name of the different variable of your environment:

		* env.workspace: Path where everything is stored (head folder). This is the path from where all the script are run
		
		* input_raster: Name of the input raster. The input raster must be projected in meters and in a projection that conserve the areas (e.g. Lambert, UTM,...) 		
		
		* raster_fnme: Name of the rasters used for the calculations (raster_fnme+"_corr", raster_fnme+"_fdir", raster_fnme+"_facc", raster_fnme+"_alti", raster_fnme+"_cost")
		
		* river_folder: Name of the generic river folders (will be increased by 1 for each next river). Generally "RProfils/R"  

 	And then, to run the script
	
	.. code-block:: python
	
		>>> execfile("1-Files-preparation.py")
 
	This produces the 5 files _corr, _fdir, _facc, _alti and _cost. They are stored in the folder "Rasters/". It also buid the structure of the workspace (see this point further)

3. In Arcgis, you have to choose the rivers you want to analyse

	- Determine the number of rivers to analyse (n)
	
	- in "RProfiles/", create a folder Ri for each river, i=[1,n] (R1/, R2/,… Rn/). R1 is already created by the 1-Files-preparation.py script. You can copy/paste it an change the number
	
	- In each "RProfiles/Ri/", you should have a "shp_file/" folder, otherwise, create it
	
	- In each "RProfiles/Ri/shp_files/", you should have 2 shapefiles, exu_Ri.shp	and src_Ri.shp, or create them (ArcCatalog: right_clic--> new_shapefile). Use ArcCatalog to rename them with the right river number (the "i")

	- In ArcMap, add the shapefiles
	
	- For each River, edit the src and exu shapefile (start_editing	--> create new feature --> use the pen to add a point). Add the source point in src_Ri.shp. Add the outlet of the river in the exu_Ri.shp

	- For each river, mesure the greatest length of the watershed (L_watershed) and write it in the file "River-caract.txt" stored in the head folder this file should be built only once:
		
		* first line =  header [river_nb,L_watershed]
		
		* then, one line per river, with 2 columns [river_nb,L_watershed]
		
		* tab separation between the columns

4. Edit the script 2-srcipt-ArcGIS.py (with a text editor)

	- Go to ####### SET your ENVIRONNEMENT  #################

	- Give the name of the different variable of your environment: It has to be consistent with the first script and your environment

		* env.workspace: Path where everything is stored (head folder). This is the path from where all the script are run

		* raster_fnme: Name of the rasters used for the calculations (raster_fnme+"_corr", raster_fnme+"_fdir", raster_fnme+"_facc", raster_fnme+"_alti", raster_fnme+"_cost")

		* river_folder: Name of the generic river folders (will be increased by 1 for each next river). Generally "RProfils/R"  

		* Rivercaract: Name of the text-file with the rivers caracteristics. first line =  header [river_nb,L_watershed], then, one line per river, with 2 columns [river_nb,L_watershed]. tab separation between the columns. Generally "River-caract.txt"

	In Arcgis, in the Python terminal, if you havn't done it (step 2), you need to go in the working folder

	.. code-block:: python

		>>> import os

		>>> os.chdir("your/path")

	And then, to run the script

	.. code-block:: python
	
		>>> execfile("2-script-ArcGIS.py")

	This script uses ArcGIS python function to build the three .dbf file for each river

5. Edit the script 3-Calcul-indices_Demoulin.py (with a text editor). It has to be consistent with the first and second scripts, as well with your environment

	- Go to ####### SET your ENVIRONNEMENT  #################
	
	- Give the name of the different variable of your environment:

		* Rivercaract: Files containing the caracteristic of the rivers you want to run: first line =  header [river_nb,L_watershed], then, one line per river, with 2 columns [river_nb,L_watershed], tab separation between the columns. Generally "River-caract.txt"
		
		* graphs_path: Folder where will be stored the graphs produced by this script. Generaly "Graphs/"                             

		* rprofiles: Folder where are stored the rivers data. Generally "RProfils/"

	- Open a terminal from where you are able to run Python scripts with the required packages
	
	- Run 

	.. code-block:: bash	
	
		python 3-Calcul-indices_Demoulin.py

	This script will produce:
		
		* Hypsometric graphs for each river stored in "Graphs/"
		
		* Regression graphs between the different rivers parameters stored in "Graphs/"
		
		* For each river, a synthesis text table "Ri-calc.txt" stored in "RProfiles/Ri"
		
		* A synthesis text table "Calc-indices-fastSR.txt" stored in the head folder
		
		* A text file containing a summary of the output ("summary.txt")


Essai_Automatic.py
------------------

There is an other script currently in development : Essai_automatic.py. 

This script will complete the script 1-Files-preparation.py. It calcules automatically:
	
	- the basins of interest based on given Strahler orders,
	
	- The outlet and source shape-files for each basin
	
	- The length of each basin
 
 and write the file River-caract.txt required for the scripts 2-script-ArcGIS.py and 3-Calcul-indices-Demoulin.py


Structure of the workspace
--------------------------

To keep a clean workspace, I built this set of scripts on a strutured project. You need to follow this structure for the database/.project:

|_RASTER

|	|_raster_fnme+"_corr"

|	|_raster_fnme+"_fdir"

|	|_raster_fnme+"_facc"

|	|_raster_fnme+"_alti"

|	|_raster_fnme+"_cost"

|_RPROFILS

|	|_Ri (i = [1:n], n = nb of rivers) 

|	|	|_TMP

|	|	|	|_temp_files_produced_by_this_script

|	|	|_SHP_FILES

|	|	|	|_exu_Ri

|	|	|	|_src_Ri

|_Shp_files

|	|_shapes_files_produced_Essai_Automatic.py

|_Graphs

|	|_Outputs_graphs

|_TMP

|	|_tmp_files

|_script_ArcGIS.py

|_Calcul-indices_Demoulin.py

|_River-caract.txt

To run a script on ArcGIS, in the ArcGIS Python console, you need to go in the working folder

.. code-block:: python

	>>> import os   #(import os module)

	>>> os.chdir("your/path")

And then, to run the script:

.. code-block:: python

	>>> execfile("script-ArcGIS.py")


Zip the results
---------------

If you have access to a Unix-like terminal (On Linux, Mac, or Win with cygwin), You can easily make a tar.gz file from your working folder by running the tarngo.sh shell script

Before to run it, edit it in a text editor, You have to change the variable and set them to your workspace
	
	- NB: Number of rivers computed

	- INPUTRASTER: Name of the original raster
	
	- DATA: Do you want input data files ? yes if yes, non if no

		DATA=yes
		
		DATA=no

To run th script, just type in your terminal :

.. code-block:: bash

	./tarngo.sh

It will produce a Calc-indices_Demoulin.tar.gz file

Good Luck

LICENCE
-------

This package is licenced with `CCby-nc <https://creativecommons.org/licenses/by-nc-sa/3.0/>`_

