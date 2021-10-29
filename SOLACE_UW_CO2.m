file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V08 SOLACE/IMOS_SOOP-CO2_GST_20201204T043818Z_VLMJ_FV01.nc';

DfCO2 = ncread(file, 'DfCO2');
DfCO2_QC = ncread(file, 'DfCO2_quality_control');
u = ncread(file,'WSPD');
lat = ncread(file, 'LATITUDE');
pressure = ncread(file, 'Press_ATM');
press_equil = ncread(file, 'Press_Equil');
time = ncread(file, 'TIME') + datetime(1950,1,1);
d = time;
doy = day(d,'dayofyear');
T = ncread(file,'TEMP');
Temp_equil = ncread(file, 'TEMP_2');
S = ncread(file,'PSAL');
sss = S;
sst = T;

% % convert from mole fraction to pCO2
% slp = press_equil./101.325;
% svp = exp(25.4543-67.4509.*(100./(sst+273.15)) - 4.8489.*log((sst+273.15)./100) - 0.000544.*sss);
% pCO2_sw_Equ = xCO2_SW.*slp.*(1-svp);
% pCO2_sw = pCO2_sw_Equ .*exp((0.0423*(sst - Temp_equil)));
% pCO2_atm = xCO2_AIR.*slp.*(1-svp);
% % Pfe2011a.pdf p.2
% % Pfeil, 2013
% %DpCO2 = pCO2_sw - pCO2_atm;


[F_CO2]=FCO2_CWE(DfCO2,T,S,u);


fig = figure()
scatter(time,F_CO2,[],lat,'filled')
c = colorbar;
c.Label.String = 'Latitude';
hold on
yline(0);
xlabel('time')
ylabel('air-sea CO2 flux SOLACE')
saveas(fig, 'SOLACE_UW_airseaFlux','png')


