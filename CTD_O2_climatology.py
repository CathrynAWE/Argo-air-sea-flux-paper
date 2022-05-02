# This scripts puts all CTD sensor data into a csv file, so I can use it later for various purposes

import pandas as pd
import os
import gsw
import netCDF4
from netCDF4 import Dataset
import numpy as np

# all CTD data is located in
CTD_folder = 'C:/Users/cawynn/cloudstor/Shared/CTD/CTD'

# compile all the necessary files
files = []
for i in os.listdir(CTD_folder):
    if i.endswith('.csv'):
        filename = [CTD_folder + '/' + i]
        files += filename

# first put all the files into one frame
CTD_data = []

for i in files:
    data = pd.read_csv(i)
    SA = gsw.SA_from_SP(data.SALINITY, data.PRESSURE, data.START_LON, data.START_LAT)
    CT = gsw.CT_from_t(SA, data.TEMPERATURE, data.PRESSURE)
    pt = gsw.pt0_from_t(SA, data.TEMPERATURE, data.PRESSURE)
    sigma0 = gsw.sigma0(SA, CT)

    # the oxygen data is presented as umol/L in the CTD csv files
    # need to convert to umol/kg
    O2_umol_kg = (data.OXYGEN / (sigma0 + 1000)) * 1000
    oxsol = gsw.O2sol_SP_pt(data.SALINITY, pt)
    doxs = O2_umol_kg / oxsol
    # add these to the dataframe before appending to the list CTD_data
    data['O2_umol_kg'] = O2_umol_kg
    data['oxsol'] = oxsol
    data['doxs'] = doxs
    CTD_data.append(data)


column_names = list(data)

# now export list to panda frame
CTD_data_ox = pd.concat(CTD_data)
# and then to csv files
CTD_data_ox.to_csv('CTD_data_ox.csv')


#### some CTD cast data is in a different format, so it's easier to just make a separate file out of those
# 2020_V08, IN2019_V02 #14, IN2018_V07 #2, IN2017_V02 #4 & #6

data =pd.DataFrame([])
CTD_nc_data=pd.DataFrame([])

for f in glob.glob('C:/Users/cawynn/cloudstor/Air sea flux manuscript/Matlab scripts/Argo-air-sea-flux-paper/CTD_nc_files/*CtdAvg.nc'):
    f = Dataset(n, mode='r', format='NETCDF4')
    data = pd.DataFrame([])

    data['OXYGEN'] = f.variables['oxygen'][0, 0, :, 0]
    data['OXYGEN_QC'] = f.variables['oxygenFlag'][0, 0, :, 0]
    data['OXYGEN_2'] = f.variables['oxygen_2'][0, 0, :, 0]
    data['OXYGEN_2_QC'] = f.variables['oxygen_2Flag'][0, 0, :, 0]

    data['SURVEY_NAME'] = [f.getncattr('Survey')] * len(data['OXYGEN'])
    data['STATION'] = [f.getncattr('Deployment')] * len(data['OXYGEN'])
    data['START_TIME'] = [f.getncattr('StartTime')] * len(data['OXYGEN'])
    data['END_TIME'] = [f.getncattr('EndTime')] * len(data['OXYGEN'])
    data['MIN_DEPTH'] = [0] * len(data['OXYGEN'])
    data['MAX_DEPTH'] = [max(f.variables['pressure'][:])] * len(data['OXYGEN'])
    data['BOTTOM_DEPTH'] = [f.getncattr('WaterDepth')] * len(data['OXYGEN'])
    data['BOTTOM_LAT'] = [f.variables['latitude'][0]] * len(data['OXYGEN'])
    data['BOTTOM_LON'] = [f.variables['longitude'][0]] * len(data['OXYGEN'])
    data['END_LAT'] = [f.variables['latitude'][0]] * len(data['OXYGEN'])
    data['END_LON'] = [f.variables['longitude'][0]] * len(data['OXYGEN'])
    data['START_LAT'] = [f.variables['latitude'][0]] * len(data['OXYGEN'])
    data['START_LON'] = [f.variables['longitude'][0]] * len(data['OXYGEN'])
    data['PRESSURE'] = f.variables['pressure'][:]

    data['SALINITY'] = f.variables['salinity'][0,0,:,0]
    data['SALINITY_QC'] = f.variables['salinityFlag'][0,0,:,0]
    data['SALINITY_2'] = f.variables['salinity_2'][0,0,:,0]
    data['SALINITY_2_QC'] = f.variables['salinity_2Flag'][0,0,:,0]

    data['TEMPERATURE'] = f.variables['temperature'][0,0,:,0]
    data['TEMPERATURE_QC'] = f.variables['temperatureFlag'][0,0,:,0]
    data['TEMPERATURE_2'] = f.variables['temperature_2'][0,0,:,0]
    data['TEMPERATURE_2_QC'] = f.variables['temperature_2Flag'][0,0,:,0]

    CTD_nc_data = pd.concat([CTD_nc_data, data], axis=0)

# and then to csv files
CTD_nc_data.to_csv('CTD_nc_data_ox.csv')