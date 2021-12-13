% air sea flux plots

% load the data
load('SOTS2020_data.mat');
load('SOTS2021_data.mat');
load('SOLACE_data.mat');
load('TEMPO_data.mat');
load('mooring_data.mat');
load('float_Data.mat');

figure()
subplot(2,1,1)
title('SOFS 2020/2021 air sea flux data')
yyaxis left
plot(SOTS2020.time(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),SOTS2020.lat(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),'-')
hold on
plot(SOTS2021.time(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),SOTS2021.lat(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),'--')
plot(SOLACE.time(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),SOLACE.lat(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),'+')
plot(TEMPO.time(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),TEMPO.lat(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),'*')
plot(float_data.time,float_data.lat)
ylabel('Latitude')

yyaxis right
plot(SOTS2020.time(SOTS2020.lon>140 & SOTS2020.lon<146),SOTS2020.lon(SOTS2020.lon>140 & SOTS2020.lon<146),'-')
plot(SOTS2021.time(SOTS2021.lon>140 & SOTS2021.lon<146),SOTS2021.lon(SOTS2021.lon>140 & SOTS2021.lon<146),'--')
plot(SOLACE.time(SOLACE.lon>140 & SOLACE.lon<146),SOLACE.lon(SOLACE.lon>140 & SOLACE.lon<146),'+')
plot(TEMPO.time(TEMPO.lon>140 & TEMPO.lon<146),TEMPO.lon(TEMPO.lon>140 & TEMPO.lon<146),'*')
plot(float_data.time,float_data.lon)
ylabel('Longitude')
legend('SOTS2020', 'SOTS2021','SOLACE','TEMPO','float','Orientation','horizontal','Location','best')
hold off
xlabel('Time')

subplot(2,1,2)
%yyaxis left
plot(SOTS2020.time(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50), SOTS2020.FCO2(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),'ob','MarkerSize',4)
yline(0)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
hold on
plot(SOTS2021.time(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50), SOTS2021.FCO2(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),'ob','MarkerSize',4)
plot(SOLACE.time(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50), SOLACE.FCO2(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),'-m','MarkerSize',4)
plot(TEMPO.time(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50), TEMPO.FCO2(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),'+c','MarkerSize',4)
plot(float_data.time,float_data.flux,'--k','MarkerSize',4)
plot(mooring_data.dxCO2_time, mooring_data.F_CO2,'or','MarkerSize',4)
legend('SOTS2020', '0','SOTS2021', 'SOLACE', 'TEMPO', 'Float', 'SOFS','Orientation','horizontal','Location','best');
% yyaxis right
% plot(datetime(mooring_data_NCP.time,'ConvertFrom','datenum'), mooring_data_NCP.ncp_C_mgm2hr,'-r','MarkerSize',4)
% ylabel('NCP C mg m-2 hr-1')
hold off
xlabel('Time')


% just mooring and SOLACE UW and SOTS2021 voyage UW
figure()
subplot(2,1,1)
title('SOFS 2020/2021 air sea flux data')
yyaxis left
plot(SOLACE.time,SOLACE.lat,'-')
hold on
plot(SOTS2021.time,SOTS2021.lat,'+','MarkerSize',4)
plot(mooring_data.dxCO2_time,mooring_data.wsp_lat,'o','MarkerSize',2)
ylabel('Latitude')

yyaxis right
plot(SOLACE.time,SOLACE.lon,'-')
plot(SOTS2021.time,SOTS2021.lon,'+','MarkerSize',2)
plot(mooring_data.dxCO2_time,mooring_data.wsp_lon,'o','MarkerSize',2)
ylabel('Longitude')
legend('SOLACE','SOTS2021','SOFS','Orientation','horizontal','Location','best')
xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
hold off
xlabel('Time')

subplot(2,1,2)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
plot(SOLACE.time, SOLACE.FCO2,'-m','MarkerSize',2)
hold on
yline(0)
plot(mooring_data.dxCO2_time, mooring_data.F_CO2,'or','MarkerSize',4)
plot(SOTS2021.time, SOTS2021.FCO2,'-c','MarkerSize',4)
legend('SOLACE','0-line','SOFS','SOTS2021','Orientation','horizontal','Location','best');
hold off
xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
xlabel('Time')



% just mooring and TEMPO UW
figure()
subplot(2,1,1)
title('SOFS 2020/2021 air sea flux data')
yyaxis left
plot(TEMPO.time,TEMPO.lat,'+')
ylabel('Latitude')
yyaxis right
plot(TEMPO.time,TEMPO.lon,'+')
ylabel('Longitude')
legend('TEMPO','Orientation','horizontal','Location','best')
hold off
xlabel('Time')

subplot(2,1,2)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
plot(TEMPO.time, TEMPO.FCO2,'-m','MarkerSize',4)
hold on
yline(0)
plot(mooring_data.dxCO2_time, mooring_data.F_CO2,'or','MarkerSize',4)
legend('TEMPO','0-line','SOFS','Orientation','horizontal','Location','best');
hold off
xlim([datetime('2021-01-29','InputFormat','yyyy-MM-dd') datetime('2021-03-26','InputFormat','yyyy-MM-dd')]);
ylim([-70 20])
xlabel('Time')
