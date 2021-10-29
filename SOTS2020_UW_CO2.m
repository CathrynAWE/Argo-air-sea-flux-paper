file = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V09 SOTS/IMOS_SOOP-CO2_GST_20200827T030011Z_VLMJ_FV01.nc';
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



[F_CO2]=FCO2_CWE(DfCO2,T,S,u);


fig = figure()
scatter(time,F_CO2,[],lat,'filled')
c = colorbar;
c.Label.String = 'Latitude';
hold on
yline(0);
xlabel('time')
ylabel('air-sea CO2 flux SOTS 2020')
ylim([-50 30]);
saveas(fig, 'SOTS2020_UW_airseaFlux','png')