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

% % read the ERA5 data for the area of the float
% file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript';
% fileERA = [file_path '\voyages_adaptor.mars.internal-1639716493.7651997-10938-14-37681699-ea52-48c7-a23c-fae4cd4304a3.nc'];
% 
% u_10 = ncread(fileERA,'u10');
% v_10 = ncread(fileERA, 'v10');
% t_2 = ncread(fileERA, 't2m')-273.16; % Kelvin converted to Celsius
% msp = ncread(fileERA, 'msl');
% EraTime = hours(ncread(fileERA, 'time'))+ datetime(1900,1,1);
% EraLat = ncread(fileERA, 'latitude');
% EraLon = ncread(fileERA, 'longitude');
% 
% wsp_TEMPO=[];
% for i = 1:length(time_TEMPO)
%     
%     % find the ERA file index that is closest in time to the UW time points
%     idx_T = knnsearch(datenum(EraTime(:)),datenum(time_TEMPO(i)));
%     idx_lat = knnsearch(EraLat(:),lat_TEMPO(i));
%     idx_lon = knnsearch(EraLon(:),lon_TEMPO(i));
% 
%     wsp_TEMPO(i) = sqrt(u_10(idx_lon,idx_lat,1,idx_T)^2 + v_10(idx_lon,idx_lat,1,idx_T)^2);
% %     t_2_f = t_2(idx_lon,idx_lat,1,idx_T);
% %     msp_f = msp(idx_lon,idx_lat,1,idx_T);
%       
% end
% 
% save('TEMPO_ERA5_Winds.mat','wsp_TEMPO','-v7.3');

load('TEMPO_ERA5_Winds.mat');


[F_CO2_TEMPO]=FCO2_CWE(DfCO2_TEMPO,T_TEMPO,S_TEMPO,wsp_TEMPO');


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