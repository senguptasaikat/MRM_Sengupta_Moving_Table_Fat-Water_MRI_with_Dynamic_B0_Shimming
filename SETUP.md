MRM Sengupta CMT Golden Angle
=============================

Requirements:
-------------
* Python - for image reconstruction
* Julia - for image reconstruction
* MATLAB - for fat/water separation (3-point, QPBO), dynamic shiimming and figure creation
* fw_i3cm0i_3point_berglund.m and associated files from ISMRM Fat Water Toolbox available at [http://ismrm.org/workshops/FatWater12/data.htm](http://ismrm.org/workshops/FatWater12/data.htm) fw_i3cm0i_3point_berglund.m is also available as FattyRiot_fw_i3cm0i_3point_berglund.m within [https://github.com/welcheb/FattyRiot](https://github.com/welcheb/FattyRiot)
* fw_i3cm1i_3pluspoint_berglund_QPBO.m and assocaited files available from [https://github.com/welcheb/fw_i3cm1i_3pluspoint_berglund_QPBO](https://github.com/welcheb/fw_i3cm1i_3pluspoint_berglund_QPBO)

Tested Configuration:
---------------------
* Mac OS X 10.10.1 (Yosemite)
* Python 2.7.8 Anaconda 2.0.1 (x86_64)
* Julia 0.3.3
	- requires NFFT and Calculus packages, i.e. Pkg.add("NFFT") and Pkg.add("Calculus")
* MATLAB R0214b

Installation Options:
---------------------
* Click the `Download ZIP` button on the lower right hand side of the repository to download the code to your local machine
* OR clone the git repository to your local machine
* OR fork to your own GitHub repository and then clone the git repository to your local machine

Usage Steps:
------------
* #### Data download
	- Run `./code/download_data_all.m` or `./code/download_data_all.py` to download the raw CMT *k*-space data files in MATLAB .mat V5 format to the `../data_input/` folder.

* #### Reconstruction, fat/water separation and dynamic shimming performed at the MRI scanner
 	- Edit variable `fileno` varibale in MATLAB m-file `./code/JULIAImageRecon_FWSeparation_DynamicShimming.m`, e.g. for subject 2 non-shimmed `fileno = 'Sub2_NS';`
	- Run MATLAB m-file `./code/JULIAImageRecon_FWSeparation_DynamicShimming.m` cell-by-cell or complete script at once to generate reconstructions of each echo, fat/water separations including off-resonance fieldmap. Output data files are written to the `./data_output/` folder.
	- Dynamic shimming is performed using a GUI tool that automatically opens follwing fat/water separation. Click `Start` button to load the fieldmaps. Click `Shim!` button to calculate shim values. Five dynamic shimming text files are written to `./data_output/` including `Shim_Table_F0.txt`, `Shim_Table_X.txt`, `Shim_Table_Y.txt`, `Shim_Table_Z.txt` and `Shim_Table_MAXVALS.txt`. Note: repeated executions of the dynamic shimming GUI overwrite these 5 files.

* #### Reconstruction and fat/water separation performed away from the MRI scanner
	- Run `./code/Python Recon/CMT_Recon.py` to generate 16 reconstructions total - 4 echoes for non-shimmed (NS) and dynamic shimmed (DS) for each subject. Output data files are written to the `./data_output/` folder.
	- Run `./code/FWSeparation_QPBO_PYTHON.m` after editing `fileno` varibale in MATLAB m-file, e.g. `fileno = 'Sub1_DS';`. Output data files are written to the `./data_output/` folder.

* #### Figure creation
	- Once results are available in `./data_output/`, figures appearing in the published paper can be recreated using any one of the scripts found in `./code` with the prefix `CMT_FWMRI_DYNSHIM_Figure`, e.g. `./code/CMT_FWMRI_DYNSHIM_Figure_1.m`
	 

Notes:
------
* The scripts `./code/download_data_verify_all.py` and `./code/download_data_verify_all.m` verify the MD5 hash values for all `./data_input/` files.
* Files in `./data_input/`, `./data_output/`, and `./figures/` (other than `README.md`) are ignored when performing `git commit`.