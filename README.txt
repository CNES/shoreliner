Modular Tool for Waterline Extraction (Graffin M., Touzé T., Bergsma E.)

Inputs : contains the config file, open and edit it whenever you want to extract waterlines on a new site or using a new method.*
### launcher.py needs files in pickle format, with a specific format (see example) to launch.

Codes : contains scripts and functions to extract the waterline.

Projects : contains all your current or past waterline extraction projects, with codes, config files and outputs (will be created the first time you run the laucnher).

### You need a certified Google Earth Engine account to use this tool. Check on https://earthengine.google.com/ to create an account, and for details.
### Create a Google Earth Engine Project (required) : https://console.cloud.google.com/projectcreate

PREWORK : 
1 - Create a conda environment with basic packages : "conda create -n shoreliner spyder numpy scipy pandas matplotlib sympy cython" (shoreliner is a proposed name for the env, feel free to change it)
2 - Activate the environment (eg, conda activate shoreliner)
3 - Download all the packages using first "pip install gdal" then "pip install -r PATH/requirements.txt"
4 - After the account has been validated, you would need to authenticate using the command "earthengine authenticate"
  OR launch spyder (type "spyder" in the shell) and type in the terminal type "import ee" then "ee.Authenticate()"
  OR use the setupGEEapi.ipynb jupyter notebook (for CNES-HPC use only).

WATERLINE EXTRACTION : 
1 - Edit the config.yaml file in ShorelinerTool/Inputs to design the extraction process as you intend it to be (site, date range, methods, ...), add your Earth Engine project ID in the config (at EE_ID), the Earth Engine project ID can be found at https://console.cloud.google.com/home/dashboard.
2 - Launch launcher.py in the ShorelinerTool folder (or sbatch job_launcher.slurm on the CNES-HPC).
3 - Go on ShorelinerTool/Projects/NAMEOFYOURPROJECT folder, you will find poly_0 to poly_n (n being the number of sub ROIs) folders, for each of them you'll have the same 3 scripts in : download.py (download images from Google Earth Engine), main_process.py (extract waterline from images and save it). Launch first download.py, then main_process.py (and optionally postprocess.py).

***IF YOU NEED TO GENERATE TRANSECTS ALONG A LARGE AREA, 
YOU CAN DOWNLOAD THE GSHHG REFERENCE COASTLINE AT https://www.soest.hawaii.edu/pwessel/gshhg/,
PLACE THE SHAPEFILES IN A SPECIFIC FOLDER AND ENTER THE CORRESPONDING PATH IN THE CONFIG FILE

