file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2021_V02 SOTS/IMOS_SOOP-CO2_GST_20210414T121956Z_VLMJ_FV01.nc';
DfCO2_SOTS2021 = ncread(file, 'DfCO2');
DfCO2_QC_SOTS2021 = ncread(file, 'DfCO2_quality_control');
u_SOTS2021 = ncread(file,'WSPD'); % this is at 24.7m and needs correcting to 10m, as per Sutton, 2017
u10_SOTS2021 = u_SOTS2021/(1+(sqrt(0.0011)/0.4)*log(24.7/10));

lat_SOTS2021 = ncread(file, 'LATITUDE');
lon_SOTS2021 = ncread(file, 'LONGITUDE');
pressure_SOTS2021 = ncread(file, 'Press_ATM');
press_equil_SOTS2021 = ncread(file, 'Press_Equil');
time_SOTS2021 = ncread(file, 'TIME') + datetime(1950,1,1);
d_SOTS2021 = time_SOTS2021;
doy_SOTS2021 = day(d_SOTS2021,'dayofyear');
T_SOTS2021 = ncread(file,'TEMP');
Temp_equil_SOTS2021 = ncread(file, 'TEMP_2');
S_SOTS2021 = ncread(file,'PSAL');
sss_SOTS2021 = S_SOTS2021;
sst_SOTS2021 = T_SOTS2021;


[F_CO2_SOTS2021]=FCO2_CWE(DfCO2_SOTS2021,T_SOTS2021,S_SOTS2021,u10_SOTS2021);



SOTS2021.FCO2 = F_CO2_SOTS2021;
SOTS2021.time = time_SOTS2021;
SOTS2021.lat = lat_SOTS2021;
SOTS2021.lon = lon_SOTS2021;
SOTS2021.TEMP = T_SOTS2021;
SOTS2021.PSAL = S_SOTS2021;
SOTS2021.dfCO2 = DfCO2_SOTS2021;
SOTS2021.dfCO2_QC = DfCO2_QC_SOTS2021;

clearvars -except SOTS2021

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('SOTS2021_data.mat')

% fig = figure()
% scatter(time_SOTS2021,F_CO2_SOTS2021,[],lat_SOTS2021,'filled')
% c = colorbar;
% c.Label.String = 'Latitude';
% hold on
% yline(0);
% xlabel('time')
% ylabel('air-sea CO2 flux SOTS 2021')
% ylim([-100 20]);
% saveas(fig, 'SOTS2021_UW_airseaFlux','png')

% figure()
% geoscatter(lat_SOTS2021, lon_SOTS2021, datenum(time_SOTS2021)/10000, F_CO2_SOTS2021,'.')
% c=colorbar;
% c.Label.String = 'Air sea flux';
% title('SOTS 2021 UW data')
% caxis([-70 30])
% 
% figure()
% subplot(2,1,1)
% title('SOTS 2021 UW data')
% yyaxis left
% plot(time_SOTS2021,lat_SOTS2021,'-r')
% ylabel('Latitude')
% yyaxis right
% plot(time_SOTS2021,lon_SOTS2021,'-b')
% ylabel('Longitude')
% xlabel('Time')
% 
% subplot(2,1,2)
% yyaxis left
% plot(time_SOTS2021,F_CO2_SOTS2021,'ok','MarkerSize',2)
% ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
% yyaxis right
% plot(time_SOTS2021,T_SOTS2021,'-b')
% ylabel('SST')
% xlabel('Time')
