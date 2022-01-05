% this script calculates SOFS mooring air-sea fluxes

% calculate flux from dfCO2, SST and S in the files

% translate windspeed from whichever mooring height to 10m
% with Sutton, 2017 equation
% u10[m s-1] = u [m s-1] / (1 + (sqrt(0.0011)/0.4) * ln(Z/10))
% Z is wind anemometer height [m]
% not needed, windspeed in files is 10m converted

% patch together daily files, with hourly obs for windspeed and flags
% http://thredds.aodn.org.au/thredds/catalog/IMOS/DWM/ASFS/SOFS/Surface_fluxes/Real-time/2020_daily/catalog.html

addpath ('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper')

path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\SOFS xCO2';

cd(path)

nc_files = dir('*.nc');

mooring_data= [];
mooring_data.wsp =[];
mooring_data.wsp_qc =[];
mooring_data.wsp_time =[];
mooring_data.lat =[];
mooring_data.lon =[];

for i = 1:length(nc_files)
    
    file = [nc_files(i).folder '\' nc_files(i).name];

    wsp = ncread(file,'WSPD10M');
    wsp_qc = ncread(file,'WIND_FLAG');
    wsp_time = ncread(file, 'TIME') + datetime(1950,1,1);
    lat = ncread(file, 'LATITUDE');
    lon = ncread(file, 'LONGITUDE');
    
    mooring_data.wsp = [mooring_data.wsp; wsp(:)];
    mooring_data.wsp_qc = [mooring_data.wsp_qc; wsp_qc(:)];
    mooring_data.wsp_time = [mooring_data.wsp_time; wsp_time(:)];
    mooring_data.lat = [mooring_data.lat; lat(:)];
    mooring_data.lon = [mooring_data.lon; lon(:)];
   
end

% download data from
% https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0118546/
% stitch CSV files together,

% in the interim I will use Pete's gridded data file

path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\SOFS xCO2\pCO2';
cd(path)

% pCO2_file = dir('*pCO2*.nc');
pCO2_file = dir('mooring*merge*.nc');


pCO2_filename = [pCO2_file.folder '\' pCO2_file.name];

xCO2_SW = ncread(pCO2_filename,'xCO2_SW');

xCO2_AIR = ncread(pCO2_filename,'xCO2_AIR');

dxCO2 = xCO2_SW - xCO2_AIR;

xCO2_Time = ncread(pCO2_filename, 'TIME') + datetime(1950,1,1);

xCO2_PSAL = ncread(pCO2_filename,'PSAL');

xCO2_SST = ncread(pCO2_filename,'TEMP');

pres = ncread(pCO2_filename, 'pressure');

%%%%%%%%%%%%%%%%%%%
temp_air = ncread(pCO2_filename, 'AIRT');

% load the mooring pCO2 data into the mooring_data array, starting from the
% first data point in 2020

mooring_data.dxCO2 = dxCO2(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_sw = xCO2_SW(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_air = xCO2_AIR(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_time = xCO2_Time(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_PSAL = xCO2_PSAL(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_SST = xCO2_SST(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_pres = pres(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_temp_air = temp_air(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));


% convert mole fraction in dry air into pCO2
% for seawater xCO2
% pCO2_sw = xCO2_SW (umol/mol) * (p - pH2O) % result in atm
% pH20_sw = exp(24.4543 - (6745.09/((T_C)+273.16)) - 4.8489 *
% log(((T_C)+273.16)/100) - 0.000544*PSAL) % result in atm, T converted to
% K

% need to translate p into atm from kPa
mooring_data.pH20_sw = exp(24.4543 - (6745.09./(mooring_data.xCO2_SST+273.16)) - ...
    4.8489*log((mooring_data.xCO2_SST+273.16)/100) - 0.000544.*mooring_data.xCO2_PSAL);

mooring_data.pCO2_sw = mooring_data.xCO2_sw .* ((mooring_data.xCO2_pres.*0.009869233) - mooring_data.pH20_sw);


% for air xCO2
% pCO2_air = xCO2_Air (umol/mol) * (p - pH2O) % result in atm
% pH20_air = 6.11 * 10^((7.5*T)/(237.3+T)) % result in mbar or hPa, T in C
% pH20_air converted to kPa = pH20/10  

mooring_data.pH20_air = (6.11 * 10.^((7.6 .* mooring_data.xCO2_temp_air)./(273.3+mooring_data.xCO2_temp_air)))/10;
% result in kPa

% need to translate delta p into atm from kPa by multiplying with*0.009869233
mooring_data.pCO2_air = mooring_data.xCO2_air .* ((mooring_data.xCO2_pres - ...
    mooring_data.pH20_air)*0.009869233); % result in atm


mooring_data.dpCO2 = mooring_data.pCO2_sw - mooring_data.pCO2_air; 

% now limit the windspeed readings to the number of xCO2 readings

l=[];
for i = 1:length(mooring_data.xCO2_time)
    l(i) = knnsearch(datenum(mooring_data.wsp_time(:)),datenum(mooring_data.xCO2_time(i)));
end

mooring_data.wsp_pCO2 = mooring_data.wsp(l);
mooring_data.wsp_pCO2_time = mooring_data.wsp_time(l);
mooring_data.wsp_lat = mooring_data.lat(l);
mooring_data.wsp_lon = mooring_data.lon(l);


mooring_data.flux_xCO2=FCO2_CWE(mooring_data.dxCO2,mooring_data.xCO2_SST,mooring_data.xCO2_PSAL,mooring_data.wsp_pCO2);
mooring_data.flux_pCO2=FCO2_CWE(mooring_data.dpCO2,mooring_data.xCO2_SST,mooring_data.xCO2_PSAL,mooring_data.wsp_pCO2);



clearvars -except mooring_data

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('mooring_data.mat')
