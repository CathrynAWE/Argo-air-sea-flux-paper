% play with windspeed averages, like weighted averages, daily averages, etc

file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V08 SOLACE/IMOS_SOOP-CO2_GST_20201204T043818Z_VLMJ_FV01.nc';

DfCO2_SOLACE = ncread(file, 'DfCO2');
DfCO2_QC_SOLACE = ncread(file, 'DfCO2_quality_control');
u_SOLACE = ncread(file,'WSPD'); % this is at 24.7m and needs correcting
% to 10m, as per Sutton, 2017, but also, this is from the ship, so I will
% use ERA5 winds instead
u10_ship_SOLACE = u_SOLACE/(1+(sqrt(0.0011)/0.4)*log(24.7/10));
lat_SOLACE = ncread(file, 'LATITUDE');
lon_SOLACE = ncread(file, 'LONGITUDE');
pressure_SOLACE = ncread(file, 'Press_ATM');
press_equil_SOLACE = ncread(file, 'Press_Equil');
time_SOLACE = ncread(file, 'TIME') + datetime(1950,1,1);
d_SOLACE = time_SOLACE;
doy_SOLACE = day(d_SOLACE,'dayofyear');
T_SOLACE = ncread(file,'TEMP'); % degrees C
Temp_equil_SOLACE = ncread(file, 'TEMP_2');
S_SOLACE = ncread(file,'PSAL');
sss_SOLACE = S_SOLACE;
sst_SOLACE = T_SOLACE;

% read the ERA5 data for the area of the float
% file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript';
% fileERA = [file_path '\adaptor.mars.internal-1635986988.5401623-6576-5-bb73848d-a029-46b4-abef-e1d94760fbeb.nc'];
% 
% u_10 = ncread(fileERA,'u10');
% v_10 = ncread(fileERA, 'v10');
% t_2 = ncread(fileERA, 't2m')-273.16; % Kelvin converted to Celsius
% msp = ncread(fileERA, 'msl');
% EraTime = hours(ncread(fileERA, 'time'))+ datetime(1900,1,1);
% EraLat = ncread(fileERA, 'latitude');
% EraLon = ncread(fileERA, 'longitude');
% 
% wsp_SOLACE=[];
% for i = 1:length(time_SOLACE)
%     
%     % find the ERA file index that is closest in time to the UW time points
%     idx_T = knnsearch(datenum(EraTime(:)),datenum(time_SOLACE(i)));
%     idx_lat = knnsearch(EraLat(:),lat_SOLACE(i));
%     idx_lon = knnsearch(EraLon(:),lon_SOLACE(i));
% 
%     wsp_SOLACE(i) = sqrt(u_10(idx_lon,idx_lat,1,idx_T)^2 + v_10(idx_lon,idx_lat,1,idx_T)^2);
% %     t_2_f = t_2(idx_lon,idx_lat,1,idx_T);
% %     msp_f = msp(idx_lon,idx_lat,1,idx_T);
%       
% end

% save('SOLACE_ERA5_Winds.mat','wsp_SOLACE','-v7.3');

load('SOLACE_ERA5_winds.mat');

[F_CO2_SOLACE]=FCO2_CWE(DfCO2_SOLACE,T_SOLACE,S_SOLACE,wsp_SOLACE');

SOLACE.FCO2 = F_CO2_SOLACE;
SOLACE.time = time_SOLACE;
SOLACE.lat = lat_SOLACE;
SOLACE.lon = lon_SOLACE;
SOLACE.TEMP = T_SOLACE;
SOLACE.PSAL = S_SOLACE;
SOLACE.dfCO2 = DfCO2_SOLACE;
SOLACE.dfCO2_QC = DfCO2_QC_SOLACE;

clearvars -except SOLACE

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('SOLACE_data.mat')

% fig = figure();
% scatter(time_SOLACE,F_CO2_SOLACE,[],lat_SOLACE,'filled')
% c = colorbar;
% c.Label.String = 'Latitude';
% hold on
% yline(0);
% xlabel('time')
% ylabel('air-sea CO2 flux SOLACE')
% saveas(fig, 'SOLACE_UW_airseaFlux','png')

% figure()
% geoscatter(lat_SOLACE, lon_SOLACE, datenum(time_SOLACE)/6000, F_CO2_SOLACE,'.')
% c=colorbar;
% c.Label.String = 'Air sea flux';
% title('SOLACE UW data')
% caxis([-70 30])

% windspeeds from ship
%[F_CO2_ship_SOLACE]=FCO2_CWE(DfCO2_SOLACE,T_SOLACE,S_SOLACE,u10_ship_SOLACE);
% 
% figure()
% subplot(2,1,1)
% title('SOLACE UW data')
% yyaxis left
% plot(time_SOLACE,lat_SOLACE,'-b')
% ylabel('Latitude')
% yyaxis right
% plot(time_SOLACE,lon_SOLACE,'-r')
% ylabel('Longitude')
% xlabel('Time')
% 
% subplot(2,1,2)
% yyaxis left
% % plot(time_SOLACE,F_CO2_ship_SOLACE,'ob','MarkerSize',2)
% % ylabel('Air-sea CO2 flux mmol m^-2 d^-1 - ships winds')
% plot(time_SOLACE,T_SOLACE,'ob','MarkerSize',2)
% ylabel('SST UW - C')
% yyaxis right
% plot(time_SOLACE,F_CO2_SOLACE,'-r')
% ylabel('Air-sea CO2 flux mmol m^-2 d^-1 - ERA5 winds')
% xlabel('Time')
