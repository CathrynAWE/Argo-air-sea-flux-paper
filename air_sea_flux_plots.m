% air sea flux plots
% close all

% load the data
load('SOTS2020_data.mat');
load('SOTS2021_data.mat');
load('SOLACE_data.mat');
load('TEMPO_data.mat');
load('mooring_data.mat');
load('SOTS_float_Data.mat');
load('S55_float_Data.mat');

% mooring vs ship and SOTS float within 1 degree lat and lon
figure()
plot(SOTS2020.time(abs(SOTS2020.lat-(-47))<=1 & abs(SOTS2020.lon-142) <=1),SOTS2020.FCO2(abs(SOTS2020.lat-(-47))<=1 & abs(SOTS2020.lon-142) <=1),'-')
hold on
yline(0)
plot(SOTS2021.time(abs(SOTS2021.lat-(-47))<=1 & abs(SOTS2021.lon-142) <=1),SOTS2021.FCO2(abs(SOTS2021.lat-(-47))<=1 & abs(SOTS2021.lon-142) <=1),'.')
plot(SOLACE.time(abs(SOLACE.lat-(-47))<=1 & abs(SOLACE.lon-142) <=1),SOLACE.FCO2(abs(SOLACE.lat-(-47))<=1 & abs(SOLACE.lon-142) <=1),'*')
plot(TEMPO.time(abs(TEMPO.lat-(-47))<=1 & abs(TEMPO.lon-142) <=1),TEMPO.FCO2(abs(TEMPO.lat-(-47))<=1 & abs(TEMPO.lon-142) <=1),'+')
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'or','MarkerSize',4)
plot(SOTS_float_data.time(abs(SOTS_float_data.lat-(-47))<=1 & abs(SOTS_float_data.lon-142) <=1),SOTS_float_data.flux_L_S(abs(SOTS_float_data.lat-(-47))<=1 & abs(SOTS_float_data.lon-142) <=1),'ok','MarkerSize',4)
legend('SOTS2020','0line','SOTS2021','SOLACE','TEMPO','SOFS','SOTS float LS','Orientation','horizontal','Location','bestoutside')
title('mooring vs ship and float flux data 47+-1 lat and 142+-1 lon')
hold off
xlabel('Time')
ylabel('air sea flux mmol m^-2_ d^-^1')

%%% SOTS float
SOLACE_SOTS_lat_msk = SOLACE.lat>min(SOTS_float_data.lat) & SOLACE.lat<max(SOTS_float_data.lat);
SOLACE_SOTS_lon_msk = SOLACE.lon>min(SOTS_float_data.lon) & SOLACE.lon<max(SOTS_float_data.lon);
SOLACE_SOTS_lat_lon_msk = SOLACE_SOTS_lat_msk+SOLACE_SOTS_lon_msk;
TEMPO_SOTS_lat_msk = TEMPO.lat>min(SOTS_float_data.lat) & TEMPO.lat<max(SOTS_float_data.lat);
TEMPO_SOTS_lon_msk = TEMPO.lon>min(SOTS_float_data.lon) & TEMPO.lon<max(SOTS_float_data.lon);
TEMPO_SOTS_lat_lon_msk = TEMPO_SOTS_lat_msk+TEMPO_SOTS_lon_msk;
SOTS2020_lat_msk = SOTS2020.lat>min(SOTS_float_data.lat) & SOTS2020.lat<max(SOTS_float_data.lat);
SOTS2020_lon_msk = SOTS2020.lon>min(SOTS_float_data.lon) & SOTS2020.lon<max(SOTS_float_data.lon);
SOTS2020_lat_lon_msk = SOTS2020_lat_msk+SOTS2020_lon_msk;
SOTS2021_lat_msk = SOTS2021.lat>min(SOTS_float_data.lat) & SOTS2021.lat<max(SOTS_float_data.lat);
SOTS2021_lon_msk = SOTS2021.lon>min(SOTS_float_data.lon) & SOTS2021.lon<max(SOTS_float_data.lon);
SOTS2021_lat_lon_msk = SOTS2021_lat_msk+SOTS2021_lon_msk;

% mooring vs ship and SOTS float within min and max of float lat and lon
figure()
plot(SOTS2020.time(SOTS2020_lat_lon_msk==2),SOTS2020.FCO2(SOTS2020_lat_lon_msk==2),'-')
hold on
yline(0)
plot(SOTS2021.time(SOTS2021_lat_lon_msk==2),SOTS2021.FCO2(SOTS2021_lat_lon_msk==2),'.')
plot(SOLACE.time(SOLACE_SOTS_lat_lon_msk==2),SOLACE.FCO2(SOLACE_SOTS_lat_lon_msk==2),'*')
plot(TEMPO.time(TEMPO_SOTS_lat_lon_msk==2),TEMPO.FCO2(TEMPO_SOTS_lat_lon_msk==2),'+')
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'^r','MarkerSize',2)
plot(SOTS_float_data.time,SOTS_float_data.flux_L_D,'ok','MarkerSize',4)
plot(SOTS_float_data.time,SOTS_float_data.flux_L_S,'ob','MarkerSize',4)
plot(SOTS_float_data.time,SOTS_float_data.flux_W_D,'om','MarkerSize',4)
plot(SOTS_float_data.time,SOTS_float_data.flux_L_D_ES,'.c','MarkerSize',11)
plot(SOTS_float_data.time,SOTS_float_data.flux_L_S_ES,'.r','MarkerSize',11)
plot(SOTS_float_data.time,SOTS_float_data.flux_W_D_ES,'.g','MarkerSize',11)
ylabel('air sea flux mmol m^-2_ d^-^1')

legend('SOTS2020','0line','SOTS2021','SOLACE','TEMPO','SOFS','SOTS float LD LIAR',...
    'SOTS float LS LIAR','SOTS float WD LIAR','SOTS float LD ES','SOTS float LS ES',...
    'SOTS float WD ES','Orientation','vertical','Location','westoutside')
title('mooring vs ship and float flux data within float max and min lat and lon')
xlim([datetime('01-08-2020','InputFormat','dd-MM-yyyy') datetime('01-09-2021','InputFormat','dd-MM-yyyy')])
hold off



%%% 55S float
load('SOLACE_data.mat')
load('TEMPO_data.mat')
SOLACE_55_lat_msk = SOLACE.lat>min(S55_float_data.lat) & SOLACE.lat<max(S55_float_data.lat);
SOLACE_55_lon_msk = SOLACE.lon>min(S55_float_data.lon) & SOLACE.lon<max(S55_float_data.lon);
SOLACE_55_lat_lon_msk = SOLACE_55_lat_msk+SOLACE_55_lon_msk;
TEMPO_55_lat_msk = TEMPO.lat>min(S55_float_data.lat) & TEMPO.lat<max(S55_float_data.lat);
TEMPO_55_lon_msk = TEMPO.lon>min(S55_float_data.lon) & TEMPO.lon<max(S55_float_data.lon);
TEMPO_55_lat_lon_msk = TEMPO_55_lat_msk+TEMPO_55_lon_msk;

figure()
subplot(2,1,1)
% title('ship vs 55S float flux data within lat and lon max and min of float')
yyaxis left
plot(SOLACE.time(SOLACE_55_lat_msk),SOLACE.lat(SOLACE_55_lat_msk),'*r','MarkerSize',2)
hold on
plot(TEMPO.time(TEMPO_55_lat_msk),TEMPO.lat(TEMPO_55_lat_msk),'+r','MarkerSize',2)
plot(S55_float_data.time,S55_float_data.lat,'or','MarkerSize',4)
% legend('SOLACE','TEMPO','55S float','Orientation','horizontal','Location','bestoutside')
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('13-01-2021','InputFormat','dd-MM-yyyy')]) 


yyaxis right
plot(SOLACE.time(SOLACE_55_lon_msk),SOLACE.lon(SOLACE_55_lon_msk),'*b','MarkerSize',2)
hold on
plot(TEMPO.time(TEMPO_55_lon_msk),TEMPO.lon(TEMPO_55_lon_msk),'+b','MarkerSize',2)
plot(S55_float_data.time,S55_float_data.lon,'ob','MarkerSize',4)
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('13-01-2021','InputFormat','dd-MM-yyyy')]) 

subplot(2,1,2)
plot(SOLACE.time(SOLACE_55_lat_lon_msk==2),SOLACE.FCO2(SOLACE_55_lat_lon_msk==2),'*b','MarkerSize',2)
hold on
plot(TEMPO.time(TEMPO_55_lat_lon_msk==2),TEMPO.FCO2(TEMPO_55_lat_lon_msk==2),'+b','MarkerSize',2)
plot(S55_float_data.time,S55_float_data.flux_L_D,'+r','MarkerSize',6)
plot(S55_float_data.time,S55_float_data.flux_L_S,'*g','MarkerSize',6)
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('13-01-2021','InputFormat','dd-MM-yyyy')]) 
ylabel('air sea flux mmol m^-2_ d^-^1')
% legend('SOLACE','TEMPO','55S float LD','55S float LS','Orientation','horizontal','Location','bestoutside')


% mooring vs ship UW data
figure()
subplot(2,1,1)
title('SOFS 2020/2021 air sea flux data within 45-50 lat and 140-146 lon')
yyaxis left
plot(SOTS2020.time(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),SOTS2020.lat(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),'-')
hold on
plot(SOTS2021.time(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),SOTS2021.lat(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),'^')
plot(SOLACE.time(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),SOLACE.lat(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),'+')
plot(TEMPO.time(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),TEMPO.lat(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),'*')
plot(SOTS_float_data.time,SOTS_float_data.lat,'o')
ylabel('Latitude')

yyaxis right
plot(SOTS2020.time(SOTS2020.lon>140 & SOTS2020.lon<146),SOTS2020.lon(SOTS2020.lon>140 & SOTS2020.lon<146),'-')
plot(SOTS2021.time(SOTS2021.lon>140 & SOTS2021.lon<146),SOTS2021.lon(SOTS2021.lon>140 & SOTS2021.lon<146),'^')
plot(SOLACE.time(SOLACE.lon>140 & SOLACE.lon<146),SOLACE.lon(SOLACE.lon>140 & SOLACE.lon<146),'+')
plot(TEMPO.time(TEMPO.lon>140 & TEMPO.lon<146),TEMPO.lon(TEMPO.lon>140 & TEMPO.lon<146),'*')
plot(SOTS_float_data.time,SOTS_float_data.lon,'o')
ylabel('Longitude')
legend('SOTS2020', 'SOTS2021','SOLACE','TEMPO','float','Orientation','horizontal','Location','bestoutside')
hold off
xlabel('Time')

subplot(2,1,2)
%yyaxis left
plot(SOTS2020.time(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50), SOTS2020.FCO2(abs(SOTS2020.lat)>45 & abs(SOTS2020.lat) <50),'-b','MarkerSize',4)
yline(0)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
hold on
plot(SOTS2021.time(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50), SOTS2021.FCO2(abs(SOTS2021.lat)>45 & abs(SOTS2021.lat) <50),'^b','MarkerSize',4)
plot(SOLACE.time(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50), SOLACE.FCO2(abs(SOLACE.lat)>45 & abs(SOLACE.lat) <50),'+m','MarkerSize',4)
plot(TEMPO.time(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50), TEMPO.FCO2(abs(TEMPO.lat)>45 & abs(TEMPO.lat) <50),'*c','MarkerSize',4)
% plot(SOTS_float_data.time,SOTS_float_data.flux,'ok','MarkerSize',4)
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'.r','MarkerSize',4)
legend('SOTS2020', '0','SOTS2021', 'SOLACE', 'TEMPO', 'Float', 'SOFS','Orientation','horizontal','Location','best');
hold off
xlabel('Time')


% just mooring and SOLACE UW and SOTS2021 voyage UW
figure()
subplot(2,1,1)
title('SOFS and ship air sea flux data')
yyaxis left
plot(SOLACE.time,SOLACE.lat,'-')
hold on
plot(SOTS2021.time,SOTS2021.lat,'+','MarkerSize',4)
plot(mooring_data.xCO2_time,mooring_data.xCO2_lat,'o','MarkerSize',2)
ylabel('Latitude')

yyaxis right
plot(SOLACE.time,SOLACE.lon,'-')
plot(SOTS2021.time,SOTS2021.lon,'+','MarkerSize',2)
plot(mooring_data.xCO2_time,mooring_data.xCO2_lon,'o','MarkerSize',2)
ylabel('Longitude')
legend('SOLACE','SOTS2021','SOFS','Orientation','horizontal','Location','best')
% xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
hold off
xlabel('Time')

subplot(2,1,2)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
plot(SOLACE.time, SOLACE.FCO2,'-m','MarkerSize',2)
hold on
yline(0)
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'or','MarkerSize',4)
plot(SOTS2021.time, SOTS2021.FCO2,'+c','MarkerSize',4)
legend('SOLACE','0-line','SOFS','SOTS2021','Orientation','horizontal','Location','best');
hold off
% xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
xlabel('Time')


% just mooring and SOLACE UW and SOTS2021 voyage UW with SST and PSAL
figure()
subplot(3,1,1)
title('SOFS and ship air sea flux data')
yyaxis left
plot(SOLACE.time,SOLACE.lat,'-')
hold on
plot(SOTS2021.time,SOTS2021.lat,'+','MarkerSize',4)
plot(mooring_data.xCO2_time,mooring_data.xCO2_lat,'o','MarkerSize',2)
ylabel('Latitude')

yyaxis right
plot(SOLACE.time,SOLACE.lon,'-')
plot(SOTS2021.time,SOTS2021.lon,'+','MarkerSize',2)
plot(mooring_data.xCO2_time,mooring_data.xCO2_lon,'o','MarkerSize',2)
ylabel('Longitude')
legend('SOLACE','SOTS2021','SOFS','Orientation','horizontal','Location','best')
% xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
hold off
xlabel('Time')

subplot(3,1,2)
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
plot(SOLACE.time, SOLACE.FCO2,'-','MarkerSize',2)
hold on
yline(0)
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'or','MarkerSize',4)
plot(SOTS2021.time, SOTS2021.FCO2,'+c','MarkerSize',4)
% legend('SOLACE','0-line','SOFS','SOTS2021','Orientation','horizontal','Location','bestoutside');
hold off
% xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
xlabel('Time')
ylabel('Air sea flux')

subplot(3,1,3)
yyaxis left
plot(SOLACE.time,SOLACE.PSAL,'-')
hold on
plot(SOTS2021.time,SOTS2021.PSAL,'+','MarkerSize',4)
plot(mooring_data.xCO2_time,mooring_data.xCO2_PSAL,'o','MarkerSize',2)
ylabel('PSAL')

yyaxis right
plot(SOLACE.time,SOLACE.TEMP,'-')
plot(SOTS2021.time,SOTS2021.TEMP,'+','MarkerSize',2)
plot(mooring_data.xCO2_time,mooring_data.xCO2_SST,'o','MarkerSize',2)
ylabel('TEMP')
% legend('SOLACE','SOTS2021','SOFS','Orientation','horizontal','Location','bestoutside')
% xlim([datetime('2020-12-04','InputFormat','yyyy-MM-dd') datetime('2021-01-15','InputFormat','yyyy-MM-dd')]);
legend('off')
hold off
xlabel('Time')



% just mooring and TEMPO UW
figure()
subplot(2,1,1)
title('SOFS and TEMPO UW air sea flux data')
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
plot(mooring_data.xCO2_time, mooring_data.flux_pCO2,'or','MarkerSize',4)
legend('TEMPO','0-line','SOFS','Orientation','horizontal','Location','best');
hold off
xlim([datetime('2021-01-29','InputFormat','yyyy-MM-dd') datetime('2021-03-26','InputFormat','yyyy-MM-dd')]);
ylim([-70 20])
xlabel('Time')


%%%% 2021 mooring monthly average fluxes vs float monthly average fluxes
load('SOTS_float_data.mat')
load('mooring_data.mat')

figure()
subplot(2,1,1)
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'-b')
hold on
% plot([1:12],(mooring_data.pCO2_2020_monthly_flux/1000)*365,'-b')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_mo_ave/1000)*365,'--g')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_mo_ave/1000)*365,'--c','MarkerSize',8)
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_WD_mo_ave/1000)*365,'--k')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_corr_mo_ave/1000)*365,'-or')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,'-og')
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_WD_corr_mo_ave/1000)*365,'--ok')

hold off
% legend('mooring 2021', 'mooring 2020', 'float LD 2021', 'float LS 2021',...
%     'float WD 2021','float LDcorr 2021','float LScorr 2021','float WDcorr 2021')
xlabel('Month')
ylabel('air sea flux mol m^-2 yr^-^1')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)
% legend('mooring 2021', 'float LD 2021','float LS 2021')
%title('SOFS float')

subplot(2,1,2)
yyaxis left
plot(SOTS_float_data.time,SOTS_float_data.lat-SOTS_float_data.mooring_lat,'-r')
ylim([-3.5 2.5])
hold on
yline(0)
hold off
ylabel('float lat - mooring lat')
yyaxis right
plot(SOTS_float_data.time,SOTS_float_data.lon-SOTS_float_data.mooring_lon,'-b')
ylabel('float lon - mooring lon')
ylim([-3.5 2.5])
xlim([datetime('01-01-2021', 'InputFormat','dd-MM-yyyy'), datetime('31-12-2021', 'InputFormat','dd-MM-yyyy')])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 2021 mooring monthly average fluxes vs float monthly average fluxes
load('S55_float_data.mat')

figure()
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LD_mo_ave/1000)*365,'--g')
hold on
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LS_mo_ave/1000)*365,'--c','MarkerSize',8)
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_WD_mo_ave/1000)*365,'--k')
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LD_corr_mo_ave/1000)*365,'--og')
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LS_corr_mo_ave/1000)*365,'--oc','MarkerSize',8)
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_WD_corr_mo_ave/1000)*365,'--ok')

hold off
legend('float LD 2021', 'float LS 2021',...
    'float WD 2021','float LDcorr 2021','float LScorr 2021','float WDcorr 2021')
xlabel('Month')
ylabel('air sea flux mol m^-2 yr^-^1')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)
title('55S float')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% comparing Cape Grim xCO2 dry air to mooring xCO2 dry air
figure()
title('xCO_2 in air')
plot(SOTS_float_data.time, SOTS_float_data.CG_xCO2_uatm,'*b')
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_air,'or','MarkerSize',2)
legend('Cape Grim xCO2 dry air','SOFS xCO2 dry air')
xlabel('Time')
ylabel('xCO_2 umol/mol')
title('SOTS float')

%%%% now comparing the converted pCO2 values
figure()
title('pCO_2 in air')
plot(SOTS_float_data.time, SOTS_float_data.pCO2_uatm,'*b')
hold on
plot(mooring_data.xCO2_time, mooring_data.pCO2_air,'or','MarkerSize',2)
legend('Cape Grim pCO2','SOFS pCO2')
xlabel('Time')
ylabel('pCO_2 uatm')
title('SOTS float')

%%%%%%%%%%%%%%%%%%%%%
%%%% comparison between float based pCO2_sw calculations and mooring
%%%  pCO2_sw as measured

load('mooring_data.mat')
figure()
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D,'.r','MarkerSize',6)
hold on
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S,'.g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D,'.c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_corr,'+r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_corr,'*g','MarkerSize',9)
% plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_corr,'+c','MarkerSize',6)
plot(mooring_data.xCO2_time, mooring_data.pCO2_sw,'ob','MarkerSize',2)
hold off
% title('pCO_2 sw')
xlabel('Time')
xlim([datetime('01-12-2020','inputFormat','dd-MM-yyyy') datetime('15-01-2022','inputFormat','dd-MM-yyyy')])
ylabel('pCO_2 uatm seawater')
% title('SOFS float')
legend('float pCO_2-LD-corr','float pCO_2-LS-corr',...
    'mooring pCO_2','Orientation','horizontal','Location',...
    'southoutside','FontSize',8)
% legend('float pCO_2- LD','float pH - LS','float pH - WD',...
%     'float pCO_2-LD-corr','float pCO_2-LS-corr','float pCO_2-WD-corr',...
%     'mooring pCO_2','Orientation','horizontal','Location',...
%     'southoutside','FontSize',8)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% floats + monte carlo error bars
clear all
load('SOTS_float_data.mat')
load('mooring_data.mat')
load('monte_carlo_SOTS.mat')

figure()
subplot(2,1,1)
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'--b','LineWidth',0.7)
    shadedErrorBar([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,...
        (mooring_data.pCO2_2021_SD_flux/1000)*365,'lineprops','b','patchSaturation',[0.03])
hold on
% plot([1:12],(mooring_data.pCO2_2020_monthly_flux/1000)*365,'-b')
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_mo_ave/1000)*365,'--g')
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_mo_ave/1000)*365,'--c','MarkerSize',8)
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_WD_mo_ave/1000)*365,'--k')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_corr_mo_ave/1000)*365,'-or')
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,'-og')
% plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_WD_corr_mo_ave/1000)*365,'--ok')
shadedErrorBar([1:12],(monte_carlo_SOTS.raw.flux_LS_moav/1000)*365,...
    (monte_carlo_SOTS.simulation.flux_LS_moSD/1000)*365,'lineprops','g','patchSaturation',[0.05])
shadedErrorBar([1:12],(monte_carlo_SOTS.raw.flux_LD_moav/1000)*365,...
    (monte_carlo_SOTS.simulation.flux_LS_moSD/1000)*365,'lineprops','r','patchSaturation',[0.03])
hold off
% legend('mooring 2021', 'mooring 2020', 'float LD 2021', 'float LS 2021',...
%     'float WD 2021','float LDcorr 2021','float LScorr 2021','float WDcorr 2021')
xlabel('Month')
ylabel('air sea flux mol m^-2 yr^-^1')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)
% legend('mooring 2021', 'float LD 2021','float LS 2021')
%title('SOFS float')
hold off

subplot(2,1,2)
yyaxis left
plot(SOTS_float_data.time,SOTS_float_data.lat-SOTS_float_data.mooring_lat,'-r')
ylim([-3.5 2.5])
hold on
yline(0)
hold off
ylabel('float lat - mooring lat')
yyaxis right
plot(SOTS_float_data.time,SOTS_float_data.lon-SOTS_float_data.mooring_lon,'-b')
ylabel('float lon - mooring lon')
ylim([-3.5 2.5])
xlim([datetime('01-01-2021', 'InputFormat','dd-MM-yyyy'), datetime('31-12-2021', 'InputFormat','dd-MM-yyyy')])


%%%%%%%%%%%%%%%%%%%%%%
%%%% monte carlo error envelop vs float based error envelop
load('mooring_data.mat')
load('monte_carlo_SOTS.mat')

figure()
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_corr_mo_ave/1000)*365,'-or')
hold on
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,'-og')
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'-b')
shadedErrorBar([1:12],(SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,(SOTS_float_data.flux_LS_corr_mo_SD/1000)*365,'lineprops','r')
shadedErrorBar([1:12],(SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,(SOTS_float_data.flux_LS_corr_mo_SD/1000)*365,'lineprops','g')
shadedErrorBar([1:12],(monte_carlo_SOTS.raw.flux_LS_moav/1000)*365,(monte_carlo_SOTS.simulation.flux_LS_moSD/1000)*365,'lineprops','r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% same for S55 float

clear all
load('S55_float_data.mat')
load('monte_carlo_S55.mat')

figure()
subplot(2,1,1)
% plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LD_mo_ave/1000)*365,'--g')
% plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LS_mo_ave/1000)*365,'--c','MarkerSize',8)
% plot(S55_float_data.mo_ave_month, (S55_float_data.flux_WD_mo_ave/1000)*365,'--k')
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LD_corr_mo_ave/1000)*365,'-or')
hold on
plot(S55_float_data.mo_ave_month, (S55_float_data.flux_LS_corr_mo_ave/1000)*365,'-og')
% plot(S55_float_data.mo_ave_month, (S55_float_data.flux_WD_corr_mo_ave/1000)*365,'--ok')
shadedErrorBar([1:12],(monte_carlo_S55.raw.flux_LS_moav/1000)*365,...
    (monte_carlo_S55.simulation.flux_LS_moSD/1000)*365,'lineprops','g','patchSaturation',[0.05])
shadedErrorBar([1:12],(monte_carlo_S55.raw.flux_LD_moav/1000)*365,...
    (monte_carlo_S55.simulation.flux_LS_moSD/1000)*365,'lineprops','r','patchSaturation',[0.03])
hold off
% legend('float LD 2021', 'float LS 2021','float WD 2021','float LDcorr 2021',...
% 'float LScorr 2021','float WDcorr 2021')
xlabel('Month')
ylabel('air sea flux mol m^-2 yr^-^1')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)
% legend('float LD 2021','float LS 2021','Orientation','horizontal','Location','southoutside')
% title('S55 float')
hold off

subplot(2,1,2)
yyaxis left
plot(S55_float_data.time,S55_float_data.lat,'-r')
% ylim([-3.5 2.5])
hold on
yline(0)
hold off
ylabel('float lat')
yyaxis right
plot(S55_float_data.time,S55_float_data.lon,'-b')
ylabel('float lon')
% ylim([-3.5 2.5])
xlim([datetime('01-01-2021', 'InputFormat','dd-MM-yyyy'), datetime('31-12-2021', 'InputFormat','dd-MM-yyyy')])
