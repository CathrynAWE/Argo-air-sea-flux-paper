% OSM 2022 conference figures
clear all
load('SOTS_float_data.mat')
load('CTD_data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% float pH vs SOLACE CTD cast pH
%%% with TS diagram
d =SOTS_float_data.pres <=100;
d(:,4:175)=0;

e = SOTS_float_data.pres <=1600;
e(:,4:175)=0;
figure()
subplot(1,2,1)
% plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
% hold on
% plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
% plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(SOTS_float_data.pH_LD_corr(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.r', 'MarkerSize', 6)
hold on
plot(SOTS_float_data.pH_LS_corr(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.g', 'MarkerSize', 6)

plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
% title('SOTS float pH and SOLACE CTD casts')
% legend('. red float LD corr','. greenfloat LS corr',...
%     '. blue CTD','Orientation','horizontal','Location','bestoutside')
% legend('Triangle red float LD','o green float LS','. red float LD corr','. greenfloat LS corr',...
%     '. blue CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 100])

subplot(1,2,2)
% plot(SOTS_float_data.TEMP(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.psal(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')& SOTS_float_data.pres <=100), '^r', 'MarkerSize', 3)
plot(SOTS_float_data.TEMP(d), SOTS_float_data.psal(d), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.Depth <=100), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')& CTD_data.raw_data.Depth <=100),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
% legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')


%%%%%%%%%%%%
%%%%% all float pCO2 data as calculated from the various pH corrections and
%%%%% Alkalinity estimates
load('mooring_data.mat')
load('SOTS_float_data.mat')

figure()
% plot(SOTS_float_data.time,SOTS_float_data.pCO2_uatm,'^k','MarkerSize',4)
hold on
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D,'.r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S,'.g','MarkerSize',6)
% plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D,'.c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_corr,'+r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_corr,'+g','MarkerSize',6)
% plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_corr,'+c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_ES,'*r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_ES,'*g','MarkerSize',6)
% plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_ES,'*c','MarkerSize',6)
plot(mooring_data.xCO2_time, mooring_data.pCO2_sw, 'ok','MarkerSize',2)
hold off
xlabel('Time')
xlim([datetime('01-12-2020','inputFormat','dd-MM-yyyy') datetime('15-01-2022','inputFormat','dd-MM-yyyy')])
ylabel('pCO_2 uatm')
legend('Cape Grim atm pCO_2','float pH-LD-LIAR','float pH-LS-LIAR',...
    'float pH-WD-LIAR','float pH LD corr','float pH LS corr','float pH WD corr',...
    'float pH-LD-ES', 'float pH-LS-ES','float pH-WD-ES','mooring pCO2 sw','Orientation',...
    'vertical','Location','westoutside','FontSize',8)
title('SOTS float')

%%% lat and lon differences between float and mooring
figure()
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
xlim([datetime('10-12-2020', 'InputFormat','dd-MM-yyyy'), datetime('27-01-2022', 'InputFormat','dd-MM-yyyy')])

%%% temperature and salinity traces of float and mooring
figure()
yyaxis left
plot(SOTS_float_data.time,SOTS_float_data.PSAL_20,'-r','LineWidth',1.8)
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_PSAL, '-.r')
ylabel('float salinity top 20m')

yyaxis right
plot(SOTS_float_data.time,SOTS_float_data.TEMP_20,'-b','LineWidth',1.8)
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_SST, '-.b')
legend('float', 'SOFS')
ylabel('float sw temperature top 20m')
xlim([datetime('10-12-2020', 'InputFormat','dd-MM-yyyy'), datetime('27-01-2022', 'InputFormat','dd-MM-yyyy')])


%%%%%%%%%%%%%%%%%%%%%%%%
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
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'-b','LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_mo_ave/1000)*365,'-r','LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_mo_ave/1000)*365,'-g','MarkerSize',8,'LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_ES_mo_ave/1000)*365,'-.r','LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_ES_mo_ave/1000)*365,'-.g','MarkerSize',8,'LineWidth',2)

plot(data.mo_ave_month, (data.flux_LD_corr_mo_ave/1000)*365,'-or')
plot(data.mo_ave_month, (data.flux_LS_corr_mo_ave/1000)*365,'-og')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% the same for 55S float
%%% float pH vs SOLACE CTD cast pH
%%% this time with TS diagram and lat / lon of casts and float location
load('S55_float_data.mat')
load('CTD_data_55S.mat')
figure()
subplot(1,3,1)
plot(S55_float_data.pH_LIR_Deep(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), S55_float_data.pres(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(S55_float_data.pH_LIR_Shallow(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), S55_float_data.pres(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
% plot(S55_float_data.pH_Williams_Deep(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), S55_float_data.pres(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data_55S.raw_data.pH(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')), CTD_data_55S.raw_data.Depth(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
title('55S float pH and SOLACE CTD casts')
legend('red triangle float LD','green o float LS','black + float WD','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])

subplot(1,3,2)
plot(S55_float_data.TEMP(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), S55_float_data.psal(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data_55S.raw_data.T_insitu(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')), CTD_data_55S.raw_data.Salinity_CTD(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')

subplot(1,3,3)
plot(S55_float_data.lat(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), S55_float_data.lon(:,S55_float_data.time < datetime('03-01-2021','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data_55S.raw_data.Latitude(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')), CTD_data_55S.raw_data.Longitude(CTD_data_55S.raw_data.date < datetime('02-01-2021','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Latitude')
ylabel('Longitude')
legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')


%%%%%%%%%%%%%%%

%%% 55S float
load('SOLACE_data.mat')


figure()
subplot(2,1,1)
plot(SOLACE.time,SOLACE.FCO2,'+b','MarkerSize',2)
hold on
plot(S55_float_data.time,S55_float_data.flux_L_D_corr,'+r','MarkerSize',6)
plot(S55_float_data.time,S55_float_data.flux_L_S_corr,'*g','MarkerSize',6)
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('12-01-2021','InputFormat','dd-MM-yyyy')]) 
ylabel('air sea flux mmol m^-2_ d^-^1')
% legend('SOLACE','TEMPO','55S float LD','55S float LS','Orientation','horizontal','Location','bestoutside')

subplot(2,1,2)
yyaxis left
plot(SOLACE.time,SOLACE.PSAL,'-r','MarkerSize',2)
hold on
plot(S55_float_data.time,S55_float_data.PSAL_20,'or','MarkerSize',4)
ylabel('salinity')
% legend('SOLACE','TEMPO','55S float','Orientation','horizontal','Location','bestoutside')
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('12-01-2021','InputFormat','dd-MM-yyyy')]) 

yyaxis right
plot(SOLACE.time,SOLACE.TEMP,'-b','MarkerSize',2)
hold on
plot(S55_float_data.time,S55_float_data.TEMP_20,'ob','MarkerSize',4)
ylabel('seawater temperature')
xlim([datetime('30-12-2020','InputFormat','dd-MM-yyyy') datetime('12-01-2021','InputFormat','dd-MM-yyyy')]) 



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


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% float alkalinity LIAR vs Shadwick vs SOLACE CTD cast pH
%%% with TS diagram
d =SOTS_float_data.pres <=100;
d(:,4:175)=0;

e = SOTS_float_data.pres <=1600;
e(:,4:175)=0;

figure()
subplot(1,2,1)
plot(SOTS_float_data.Alk_LIAR(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.Alk_ES(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(CTD_data.raw_data.Alkalinity(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Alk umol/kg')
ylabel('Depth dbar')
title('SOTS float')
legend('Triangle float LIAR Alk','o float ES Alk','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 100])

subplot(1,2,2)
plot(SOTS_float_data.TEMP(d), SOTS_float_data.psal(d), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.Depth <=100), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')& CTD_data.raw_data.Depth <=100),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
% legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')
