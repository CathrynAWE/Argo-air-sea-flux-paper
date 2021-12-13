file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2021_V01 TEMPO/IMOS_SOOP-CO2_GST_20210129T004303Z_VLMJ_FV01.nc';
DfCO2_TEMPO = ncread(file, 'DfCO2');
DfCO2_QC_TEMPO = ncread(file, 'DfCO2_quality_control');
u = ncread(file,'WSPD'); % this is at 24.7m and needs correcting to 10m, as per Sutton, 2017
u10_TEMPO = u/(1+(sqrt(0.0011)/0.4)*log(24.7/10));
lat_TEMPO = ncread(file, 'LATITUDE');
lon_TEMPO = ncread(file, 'LONGITUDE');
pressure = ncread(file, 'Press_ATM');
press_equil = ncread(file, 'Press_Equil');
time_TEMPO = ncread(file, 'TIME') + datetime(1950,1,1);
d = time_TEMPO;
doy = day(d,'dayofyear');
T_TEMPO = ncread(file,'TEMP');
Temp_equil = ncread(file, 'TEMP_2');
S_TEMPO = ncread(file,'PSAL');
sss = S_TEMPO;
sst = T_TEMPO;

% % convert from mole fraction to pCO2
% slp = press_equil./101.325;
% svp = exp(25.4543-67.4509.*(100./(sst+273.15)) - 4.8489.*log((sst+273.15)./100) - 0.000544.*sss);
% pCO2_sw_Equ = xCO2_SW.*slp.*(1-svp);
% pCO2_sw = pCO2_sw_Equ .*exp((0.0423*(sst - Temp_equil)));
% pCO2_atm = xCO2_AIR.*slp.*(1-svp);
% % Pfe2011a.pdf p.2
% % Pfeil, 2013
% %DpCO2 = pCO2_sw - pCO2_atm;


[F_CO2_TEMPO]=FCO2_CWE(DfCO2_TEMPO,T_TEMPO,S_TEMPO,u10_TEMPO);


TEMPO.FCO2 = F_CO2_TEMPO;
TEMPO.time = time_TEMPO;
TEMPO.lat = lat_TEMPO;
TEMPO.lon = lon_TEMPO;
TEMPO.TEMP = T_TEMPO;
TEMPO.PSAL = S_TEMPO;
TEMPO.dfCO2 = DfCO2_TEMPO;
TEMPO.dfCO2_QC = DfCO2_QC_TEMPO;

clearvars -except TEMPO

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('TEMPO_data.mat')
% 
% fig = figure()
% scatter(time_TEMPO,F_CO2_TEMPO,[],lat_TEMPO,'filled')
% c = colorbar;
% c.Label.String = 'Latitude';
% hold on
% yline(0);
% xlabel('time')
% ylabel('air-sea CO2 flux TEMPO')
% ylim([-200 20]);
% saveas(fig, 'TEMPO_UW_airseaFlux','png')