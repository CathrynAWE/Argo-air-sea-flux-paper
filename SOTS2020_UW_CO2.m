file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V09 SOTS/IMOS_SOOP-CO2_GST_20200827T030011Z_VLMJ_FV01.nc';
DfCO2_SOTS2020 = ncread(file, 'DfCO2');
DfCO2_QC_SOTS2020 = ncread(file, 'DfCO2_quality_control');
u_SOTS2020 = ncread(file,'WSPD');
lat_SOTS2020 = ncread(file, 'LATITUDE');
lon_SOTS2020 = ncread(file, 'LONGITUDE');
pressure_SOTS2020 = ncread(file, 'Press_ATM');
press_equil_SOTS2020 = ncread(file, 'Press_Equil');
time_SOTS2020 = ncread(file, 'TIME') + datetime(1950,1,1);
d_SOTS2020 = time_SOTS2020;
doy_SOTS2020 = day(d_SOTS2020,'dayofyear');
T_SOTS2020 = ncread(file,'TEMP');
Temp_equil_SOTS2020 = ncread(file, 'TEMP_2');
S_SOTS2020 = ncread(file,'PSAL');
sss_SOTS2020 = S_SOTS2020;
sst_SOTS2020 = T_SOTS2020;



[F_CO2_SOTS2020]=FCO2_CWE(DfCO2_SOTS2020,T_SOTS2020,S_SOTS2020,u_SOTS2020);


% fig = figure()
% scatter(time_SOTS2020,F_CO2_SOTS2020,[],lat_SOTS2020,'filled')
% c = colorbar;
% c.Label.String = 'Latitude';
% hold on
% yline(0);
% xlabel('time')
% ylabel('air-sea CO2 flux SOTS 2020')
% ylim([-50 30]);
% saveas(fig, 'SOTS2020_UW_airseaFlux','png')

figure()
geoscatter(lat_SOTS2020, lon_SOTS2020, datenum(time_SOTS2020)/6000, F_CO2_SOTS2020,'.')
c=colorbar;
c.Label.String = 'Air sea flux';
title('SOTS 2020 UW data')
caxis([-70 30])