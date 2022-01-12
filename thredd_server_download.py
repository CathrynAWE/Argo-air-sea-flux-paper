import requests
import shutil as shutil
import os

# open text file with SBE file names downloaded manually from AODN and manually trimmed to the files we want
files = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/SOFS xCO2/' + 'IMOS_-_Deep_Water_Moorings_-_Southern_Ocean_Flux_Station_(SOFS)_-_Surface_fluxes_(real-time)_URLs.txt'

netCDFfiles=[]
with open(files, 'r') as f:
    netCDFfiles = f.read().splitlines()


# download the relevant SBE-ODO files for plotting comparison

def download_file_from_server_endpoint(server_endpoint, local_file_path):
   # Send HTTP GET request to server and attempt to receive a response
   response = requests.get(server_endpoint)
# If the HTTP GET request can be served
   if response.status_code == 200:
        # Write the file contents in the response to a file specified by local_file_path
    with open(local_file_path, 'wb') as local_file:
        for chunk in response.iter_content(chunk_size=128):
            local_file.write(chunk)


for f in netCDFfiles:
    fnsplit = f.split('/')
    download_file_from_server_endpoint(f, fnsplit[-1])

# move the files
dest = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/SOFS xCO2'
shutil.copy(f, dest)
