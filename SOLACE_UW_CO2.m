% play with windspeed averages, like weighted averages, daily averages, etc

file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V08 SOLACE/IMOS_SOOP-CO2_GST_20201204T043818Z_VLMJ_FV01.nc';

DfCO2_SOLACE = ncread(file, 'DfCO2');
DfCO2_QC_SOLACE = ncread(file, 'DfCO2_quality_control');
u_SOLACE = ncread(file,'WSPD');
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

% % convert from mole fraction to pCO2
% slp = press_equil./101.325;
% svp = exp(25.4543-67.4509.*(100./(sst+273.15)) - 4.8489.*log((sst+273.15)./100) - 0.000544.*sss);
% pCO2_sw_Equ = xCO2_SW.*slp.*(1-svp);
% pCO2_sw = pCO2_sw_Equ .*exp((0.0423*(sst - Temp_equil)));
% pCO2_atm = xCO2_AIR.*slp.*(1-svp);
% % Pfe2011a.pdf p.2
% % Pfeil, 2013
% %DpCO2 = pCO2_sw - pCO2_atm;


[F_CO2_SOLACE]=FCO2_CWE(DfCO2_SOLACE,T_SOLACE,S_SOLACE,u_SOLACE);


fig = figure();
scatter(time_SOLACE,F_CO2_SOLACE,[],lat_SOLACE,'filled')
c = colorbar;
c.Label.String = 'Latitude';
hold on
yline(0);
xlabel('time')
ylabel('air-sea CO2 flux SOLACE')
saveas(fig, 'SOLACE_UW_airseaFlux','png')

figure()
geoscatter(lat_SOLACE, lon_SOLACE, datenum(time_SOLACE)/6000, F_CO2_SOLACE,'.')
c=colorbar;
c.Label.String = 'Air sea flux';
title('SOLACE UW data')
caxis([-70 30])