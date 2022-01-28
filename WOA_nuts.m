%%%% for the CO2SYS functions I need nutrients, here I will extract them
%%%% for the float profiles from WOA 2018

float_ID = '5906623' % SOTS float
%float_ID = '5906624' %55S float

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)
% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name];

% read necessary variables from float profile
lat_float = ncread(file,'LATITUDE');
lon_float = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
month_float = month(time_float);


addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\nctoolbox-master\nctoolbox-master'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\ocean_data_tools-master\ocean_data_tools-master\ocean_data_tools'
setup_nctoolbox

variable_list = {'silicate','phosphate'};
time = '01';
xcoords = lon_float;
ycoords = lat_float;
zgrid = 0.1;

[woa] =  woa_build_profiles(variable_list,time,xcoords,ycoords,zgrid); % zgrid optional, no interpolation if unspecified

%https://www.ncei.noaa.gov/thredds-ocean/fileServer/ncei/woa/silicate/all/1.00/woa18_all_i01_01.nc