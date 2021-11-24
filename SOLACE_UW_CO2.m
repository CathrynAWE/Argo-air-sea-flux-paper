% play with windspeed averages, like weighted averages, daily averages, etc

file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V08 SOLACE/IMOS_SOOP-CO2_GST_20201204T043818Z_VLMJ_FV01.nc';

DfCO2_SOLACE = ncread(file, 'DfCO2');
DfCO2_QC_SOLACE = ncread(file, 'DfCO2_quality_control');
u_SOLACE = ncread(file,'WSPD'); % this is at 24.7m and needs correcting to 10m, as per Sutton, 2017
u10_SOLACE = u_SOLACE/(1+(sqrt(0.0011)/0.4)*log(24.7/10));
u10_T_SOLACE = u_SOLACE *(10/24.7)^.1;
lat_SOLACE = ncread(file, 'LATITUDE');
lon_SOLACE = ncread(file, 'LONGITUDE');
pressure_SOLACE = ncread(file, 'Press_ATM');
press_equil_SOLACE = ncread(file, 'Press_Equil');
time_SOLACE = ncread(file, 'TIME') + datetime(1950,1,1);
d_SOLACE = time_SOLACE;
doy_SOLACE = day(d_SOLACE,'dayofyear');
T_SOLACE = ncread(file,'TEMP');
Temp_equil_SOLACE = ncread(file, 'TEMP_2');
S_SOLACE = ncread(file,'PSAL');
sss_SOLACE = S_SOLACE;
sst_SOLACE = T_SOLACE;


[F_CO2_SOLACE]=FCO2_CWE(DfCO2_SOLACE,T_SOLACE,S_SOLACE,u10_SOLACE);


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

[F_CO2_T_SOLACE]=FCO2_CWE(DfCO2_SOLACE,T_SOLACE,S_SOLACE,u10_T_SOLACE);
figure()
subplot(2,1,1)
title('SOLACE UW data')
yyaxis left
plot(time_SOLACE,lat_SOLACE,'-b')
ylabel('Latitude')
yyaxis right
plot(time_SOLACE,lon_SOLACE,'-r')
ylabel('Longitude')
xlabel('Time')

subplot(2,1,2)
yyaxis left
plot(time_SOLACE,F_CO2_T_SOLACE,'ob','MarkerSize',2)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1 - Tyler winds')
yyaxis right
plot(time_SOLACE,F_CO2_SOLACE,'-r')
ylabel('Air-sea CO2 flux mmol m^-2 d^-1 - Sutten winds')
xlabel('Time')
