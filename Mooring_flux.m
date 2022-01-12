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
mooring_data.wsp_lat =[];
mooring_data.wsp_lon =[];

for i = 1:length(nc_files)
    
    file = [nc_files(i).folder '\' nc_files(i).name];

    wsp = ncread(file,'WSPD10M');
    wsp_qc = ncread(file,'WIND_FLAG');
    wsp_time = ncread(file, 'TIME') + datetime(1950,1,1);
    wsp_lat = ncread(file, 'LATITUDE');
    wsp_lon = ncread(file, 'LONGITUDE');
    
    mooring_data.wsp = [mooring_data.wsp; wsp(:)];
    mooring_data.wsp_qc = [mooring_data.wsp_qc; wsp_qc(:)];
    mooring_data.wsp_time = [mooring_data.wsp_time; wsp_time(:)];
    mooring_data.wsp_lat = [mooring_data.wsp_lat; wsp_lat(:)];
    mooring_data.wsp_lon = [mooring_data.wsp_lon; wsp_lon(:)];
   
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

xCO2_Skin_temp = ncread(pCO2_filename, 'SST');

pres = ncread(pCO2_filename, 'pressure');

lat = ncread(pCO2_filename, 'LATITUDE');

lon = ncread(pCO2_filename, 'LONGITUDE');


temp_air = ncread(pCO2_filename, 'AIRT');

% load the mooring pCO2 data into the mooring_data array, starting from the
% first data point in 2020

mooring_data.dxCO2 = dxCO2(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % umol/mol

mooring_data.xCO2_sw = xCO2_SW(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % umol/mol

mooring_data.xCO2_air = xCO2_AIR(xCO2_Time>= datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % umol/mol

mooring_data.xCO2_time = xCO2_Time(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_PSAL = xCO2_PSAL(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_SST = xCO2_SST(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % Celsius
mooring_data.xCO2_SKT = xCO2_Skin_temp(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % Celsius

mooring_data.xCO2_pres = pres(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy'));

mooring_data.xCO2_temp_air = temp_air(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % Celsius

mooring_data.xCO2_lat = lat(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % Celsius
mooring_data.xCO2_lon = lon(xCO2_Time >=datetime('01-01-2020','InputFormat','dd-MM-yyyy')); % Celsius


% now limit the windspeed readings to the number of xCO2 readings

l=[];
for i = 1:length(mooring_data.xCO2_time)
    l(i) = knnsearch(datenum(mooring_data.wsp_time(:)),datenum(mooring_data.xCO2_time(i)));

    mooring_data.wsp_pCO2 = mooring_data.wsp(l);
    mooring_data.wsp_pCO2_time = mooring_data.wsp_time(l);
    mooring_data.wsp_pCO2_lat = mooring_data.wsp_lat(l);
    mooring_data.wsp_pCO2_lon = mooring_data.wsp_lon(l);

end

%%% for some reason there is no windspeed mooring data between 
% 24th April 2021 and 01st June 2021

idx = find(datenum(mooring_data.wsp_pCO2_time) == datenum(datetime('23-04-2021 23:59:00','InputFormat','dd-MM-yyyy HH:mm:ss')));
idx_start = idx(1) +1;
idx = find(datenum(mooring_data.wsp_pCO2_time) == datenum(datetime('01-06-2021 02:59:00','InputFormat','dd-MM-yyyy HH:mm:ss')));
idx_end = idx(end) -1;

mooring_data.wsp_gap_lat = mooring_data.xCO2_lat(idx_start:idx_end);
mooring_data.wsp_gap_lon = mooring_data.xCO2_lon(idx_start:idx_end);
mooring_data.wsp_gap_time = mooring_data.xCO2_time(idx_start:idx_end);


% read the ERA5 data for the area of the float
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript';
% file = [file_path '\S55_adaptor.mars.internal-1639715049.643984-2437-12-3c2e7bca-d3f8-4d73-a52d-ede241f56e62.nc'];
file = [file_path '\SOTS_float_adaptor.mars.internal-1635986988.5401623-6576-5-bb73848d-a029-46b4-abef-e1d94760fbeb.nc'];


u_10 = ncread(file,'u10');
v_10 = ncread(file, 'v10');
t_2 = ncread(file, 't2m')-273.16; % Kelvin converted to Celsius
msp = ncread(file, 'msl');
EraTime = hours(ncread(file, 'time'))+ datetime(1900,1,1);
EraLat = ncread(file, 'latitude');
EraLon = ncread(file, 'longitude');

for i = 1:length(mooring_data.wsp_gap_lat)
% for i = 1   
    % find the ERA file index that is closest in time to the float time
    idx_T = knnsearch(datenum(EraTime(:)),datenum(mooring_data.wsp_gap_time(i)));
    idx_lat = knnsearch(EraLat(:),mooring_data.wsp_gap_lat(i));
    idx_lon = knnsearch(EraLon(:),mooring_data.wsp_gap_lon(i));

    mooring_data.wsp_gap(i) = sqrt(u_10(idx_lon,idx_lat,1,idx_T)^2 + v_10(idx_lon,idx_lat,1,idx_T)^2);
    mooring_data.air_temp_gap(i) = t_2(idx_lon,idx_lat,1,idx_T);
    mooring_data.wsp_gap_ERA_time(i) = EraTime(idx_T);
end

% insert this into the wsp mooring data subset for flux calculations
mooring_data.wsp_pCO2_time(idx_start:idx_end)=mooring_data.wsp_gap_ERA_time;
mooring_data.wsp_pCO2(idx_start:idx_end)=mooring_data.wsp_gap;

mooring_data.xCO2_temp_air(idx_start:idx_end)=mooring_data.air_temp_gap;



% convert mole fraction in dry air into pCO2
% for seawater xCO2
% pCO2_sw = xCO2_SW (umol/mol) * (p - pH2O) % result in uatm
% pH20_sw = exp(24.4543 - (6745.09/((T_C)+273.16)) - 4.8489 *
% log(((T_C)+273.16)/100) - 0.000544*PSAL) % result in atm, T converted to
% K

% need to translate p into atm from kPa
mooring_data.pH2O_sw = exp(24.4543 - (6745.09./(mooring_data.xCO2_SST+273.16)) - ...
    4.8489*log((mooring_data.xCO2_SST+273.16)/100) - 0.000544.*mooring_data.xCO2_PSAL);

mooring_data.pCO2_sw = mooring_data.xCO2_sw .* ((mooring_data.xCO2_pres.*0.009869233) - mooring_data.pH2O_sw);
% result in uatm

% for air xCO2
% pCO2_air = xCO2_Air (umol/mol) * (p - pH2O) % result in uatm
% pH20_air = 6.11 * 10^((7.5*T)/(237.3+T)) % result in mbar or hPa, T in C
% pH20_air converted to kPa = pH20/10  

mooring_data.pH2O_air = (6.11 * 10.^((7.6 .* mooring_data.xCO2_temp_air)./(273.3+mooring_data.xCO2_temp_air)))/10;
% result in kPa

% need to translate delta p into atm from kPa by multiplying with*0.009869233
mooring_data.pCO2_air = mooring_data.xCO2_air .* ((mooring_data.xCO2_pres - ...
    mooring_data.pH2O_air)*0.009869233); % result in uatm


mooring_data.dpCO2 = mooring_data.pCO2_sw - mooring_data.pCO2_air; 
% uatm

%%%% now calculate fluxes with pCO2 and xCO2 (for sanity check)
%%% result in mmol m-2 d-1
mooring_data.flux_xCO2=FCO2_CWE(mooring_data.dxCO2,mooring_data.xCO2_SST,mooring_data.xCO2_PSAL,mooring_data.wsp_pCO2);
mooring_data.flux_pCO2=FCO2_CWE(mooring_data.dpCO2,mooring_data.xCO2_SST,mooring_data.xCO2_PSAL,mooring_data.wsp_pCO2);
mooring_data.flux_pCO2_SKT=FCO2_CWE(mooring_data.dpCO2,mooring_data.xCO2_SKT,mooring_data.xCO2_PSAL,mooring_data.wsp_pCO2);


mooring_data.xCO2_month = month(mooring_data.xCO2_time);
mooring_data.xCO2_year = year(mooring_data.xCO2_time);
% mooring_data.xCO2_doy = day(mooring_data.xCO2_time,'dayofyear');
mooring_data.pCO2_2020_monthly_flux = accumarray(mooring_data.xCO2_month(mooring_data.xCO2_year==2020),mooring_data.flux_pCO2(mooring_data.xCO2_year==2020),[],@(x)mean(x,'omitnan'));
mooring_data.pCO2_2020_monthly_flux(mooring_data.pCO2_2020_monthly_flux==0)=NaN;

mooring_data.pCO2_2021_monthly_flux = accumarray(mooring_data.xCO2_month(mooring_data.xCO2_year==2021),mooring_data.flux_pCO2(mooring_data.xCO2_year==2021),[],@(x)mean(x,'omitnan'));
mooring_data.pCO2_2021_monthly_flux(mooring_data.pCO2_2021_monthly_flux==0)=NaN;


clearvars -except mooring_data

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('mooring_data.mat')

% visualize the data
figure()
plot(mooring_data.xCO2_time, mooring_data.xCO2_sw,'or')
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_air,'+b')
hold off
xlabel('Time')
ylabel('xCO_2')
legend('xCO_2 SW', 'xCO_2 air')


figure()
plot(mooring_data.xCO2_time, mooring_data.dpCO2, 'or')
hold on
plot(mooring_data.xCO2_time, mooring_data.dxCO2, '*b')
hold off
xlabel('Time')
ylabel('seawater - air CO_2')
legend('dpCO_2' , 'dxCO_2')

figure()
plot(mooring_data.wsp_time, mooring_data.wsp,'or')
hold on
plot(mooring_data.wsp_gap_ERA_time, mooring_data.wsp_gap,'ob')
hold off
xlabel('Time')
ylabel('wsp')
legend('mooring' , 'ERA5')

figure()
plot(mooring_data.xCO2_time, mooring_data.dpCO2, 'or')
hold on
plot(mooring_data.xCO2_time, mooring_data.dxCO2, '*b')
plot(mooring_data.xCO2_time, mooring_data.xCO2_PSAL, '+b')
plot(mooring_data.xCO2_time, mooring_data.xCO2_SST, '+k')
plot(mooring_data.xCO2_time, mooring_data.xCO2_pres, '+c')
plot(mooring_data.xCO2_time, mooring_data.xCO2_temp_air, '+g')
hold off
xlabel('Time')
legend('dpCO_2' , 'dxCO_2', 'PSAL', 'SST','air pressure','air temp')

figure()
plot(mooring_data.xCO2_time, (mooring_data.flux_xCO2/1000)*365, 'or')
hold on
plot(mooring_data.xCO2_time, (mooring_data.flux_pCO2/1000)*365, '*b')
hold off
xlabel('Time')
ylabel('Air sea flux mol m^-^2 yr^-^1')
legend('xCO_2', 'pCO_2')

figure()
yyaxis left
plot(mooring_data.xCO2_time, mooring_data.pH2O_sw,'or')
hold on
plot(mooring_data.xCO2_time, mooring_data.pH2O_air,'ob')
hold off
ylabel('pH_2O uatm')

yyaxis right
plot(mooring_data.xCO2_time, mooring_data.pCO2_sw,'*g')
hold on
plot(mooring_data.xCO2_time, mooring_data.pCO2_air,'*m')
hold off
xlabel('Time')
ylabel('pCO_2 uatm')
legend('pH_2O SW', 'pH_2O Air', 'pCO_2 SW', 'pCO_2 air')

figure()
yyaxis left
plot([1:12],mooring_data.pCO2_2021_monthly_flux,'-r')
hold on
plot([1:12],mooring_data.pCO2_2020_monthly_flux,'-b')
hold off
legend('2021','2020')
xlabel('Month')
ylabel('air sea flux mmol m^-2 d^-^1')

yyaxis right
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'*r')
hold on
plot([1:12],(mooring_data.pCO2_2020_monthly_flux/1000)*365,'*b')
hold off
ylabel('air sea flux mol m^-2 y^-^1')

xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])