%%% since the ERA5 data has a gap from Sept to Dec 2021 (at least for now),
%%% I am using NCEP reanalysis data to patch the gap
% https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis2/surface/catalog.html
% for mean sea level pressure (mslp)
% https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/surface/catalog.html
% for windspeed (from u and v), air temperature at the surface

clear all


% float_ID ='5906623'; %SOTS
float_ID = '5906624' % 55S


path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\NOAA_NCEP_Reanalysis';
cd(path)

%%% mslp, mean sea level pressure
f = dir('mslp*.nc');
file = [f.folder '\' f.name];
mslp_ncep = ncread(file,'mslp');%Pa
mslp_ncep_time = hours(ncread(file, 'time'))+ datetime(1800,1,1);
mslp_ncep_lat = ncread(file,'lat');
mslp_ncep_lon = ncread(file,'lon');

%%% windspeed, uwnd
f = dir('uwnd*.nc');
file = [f.folder '\' f.name];
uwnd_ncep = ncread(file,'uwnd');%m s-1
uwnd_ncep_time = hours(ncread(file, 'time'))+ datetime(1800,1,1);
uwnd_ncep_lat = ncread(file,'lat');
uwnd_ncep_lon = ncread(file,'lon');

%%% windspeed, vwnd
f = dir('vwnd*.nc');
file = [f.folder '\' f.name];
vwnd_ncep = ncread(file,'vwnd');%m s-1
vwnd_ncep_time = hours(ncread(file, 'time'))+ datetime(1800,1,1);
vwnd_ncep_lat = ncread(file,'lat');
vwnd_ncep_lon = ncread(file,'lon');

%%% air temperature at the sea surface
f = dir('air*.nc');
file = [f.folder '\' f.name];
airT_ncep = ncread(file,'air')-273.15;% Kelvin to Celsius
airT_ncep_time = hours(ncread(file, 'time'))+ datetime(1800,1,1);
airT_ncep_lat = ncread(file,'lat');
airT_ncep_lon = ncread(file,'lon');


%%%% now I will load the float profiles

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)
% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name];

% read necessary variables from float profile
lat_float = ncread(file,'LATITUDE');
lon_float = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);

for i = 1:size(lat_float,1) 
% find the NCEP file index that is closes in time and space to the float
    idx_nT = knnsearch(datenum(mslp_ncep_time(:)),datenum(time_float(i)));
    idx_nlat = knnsearch(mslp_ncep_lat(:),lat_float(i));
    idx_nlon = knnsearch(mslp_ncep_lon(:),lon_float(i));
    float_ncep.mslp_ncep_f(i) = mslp_ncep(idx_nlon,idx_nlat,idx_nT,1); %Pa
    float_ncep.mslp_ncep_time_f(i) = mslp_ncep_time(idx_nT);
    
    idx_nT = knnsearch(datenum(airT_ncep_time(:)),datenum(time_float(i)));
    idx_nlat = knnsearch(airT_ncep_lat(:),lat_float(i));
    idx_nlon = knnsearch(airT_ncep_lon(:),lon_float(i));
    float_ncep.airT_ncep_f(i) = airT_ncep(idx_nlon,idx_nlat,idx_nT,1); %Kelvin
    float_ncep.airT_ncep_time_f(i) = airT_ncep_time(idx_nT);
    
    idx_nT = knnsearch(datenum(uwnd_ncep_time(:)),datenum(time_float(i)));
    idx_nlat = knnsearch(uwnd_ncep_lat(:),lat_float(i));
    idx_nlon = knnsearch(uwnd_ncep_lon(:),lon_float(i));
    float_ncep.uwnd_ncep_f(i) = uwnd_ncep(idx_nlon,idx_nlat,idx_nT,1); %m s-1
    float_ncep.uwnd_ncep_time_f(i) = uwnd_ncep_time(idx_nT);
    
    idx_nT = knnsearch(datenum(vwnd_ncep_time(:)),datenum(time_float(i)));
    idx_nlat = knnsearch(vwnd_ncep_lat(:),lat_float(i));
    idx_nlon = knnsearch(vwnd_ncep_lon(:),lon_float(i));
    float_ncep.vwnd_ncep_f(i) = vwnd_ncep(idx_nlon,idx_nlat,idx_nT,1); %m s-1
    float_ncep.vwnd_ncep_time_f(i) = vwnd_ncep_time(idx_nT);

    float_ncep.wnd_ncep_f(i) = sqrt(float_ncep.uwnd_ncep_f(i)^2 + float_ncep.vwnd_ncep_f(i)^2); % m s-1 
end


path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper';
cd(path)
save ('NCEP_S55_float.mat')

