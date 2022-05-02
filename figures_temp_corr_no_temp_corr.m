load('SOTS_float_data_TEMP.mat')
load('SOTS_float_data.mat')
load('CTD_data.mat')
cols=lines;

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
plot(SOTS_float_data.pH_LD_corr(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.', 'MarkerSize', 10,'color',cols(3,:))
hold on
plot(SOTS_float_data_TEMP.pH_LD_corr(:,SOTS_float_data_TEMP.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data_TEMP.pres(:,SOTS_float_data_TEMP.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.k')%, 'MarkerSize', 10,'color',cols(3,:))

plot(SOTS_float_data.pH_LS_corr(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.g', 'MarkerSize', 10,'color',cols(5,:))
plot(SOTS_float_data_TEMP.pH_LS_corr(:,SOTS_float_data_TEMP.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data_TEMP.pres(:,SOTS_float_data_TEMP.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '.b')%, 'MarkerSize', 10,'color',cols(5,:))

plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 20)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
% title('SOTS float pH and SOLACE CTD casts')
% legend('. red float LD corr','. greenfloat LS corr',...
%     '. blue CTD','Orientation','horizontal','Location','bestoutside')
% legend('Triangle red float LD','o green float LS','. red float LD corr','. greenfloat LS corr',...
%     '. blue CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 1600])

subplot(1,2,2)
% plot(SOTS_float_data.TEMP(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.psal(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')& SOTS_float_data.pres <=100), '^r', 'MarkerSize', 3)
plot(SOTS_float_data.TEMP(e), SOTS_float_data.psal(e), '^k', 'MarkerSize', 4)
hold on
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.Depth <=1600), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')& CTD_data.raw_data.Depth <=1600),'.b', 'MarkerSize', 20)
hold off
xlabel('Temp')
ylabel('PSAL')
% legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')



%%%%%%%%%%%%%%%%%%%%%%
%%% floats + monte carlo error bars
clear all
load('SOTS_float_data.mat')
load('SOTS_float_data_TEMP.mat')
load('mooring_data.mat')
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
cols=lines;

figure()
subplot(2,1,1)
plot([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,'.b','MarkerSize',18,'LineWidth',3)
% shadedErrorBar([1:12],(mooring_data.pCO2_2021_monthly_flux/1000)*365,...
%         (mooring_data.pCO2_2021_SD_flux/1000)*365,'lineprops','-.b','patchSaturation',[0.03])
hold on
plot(SOTS_float_data_TEMP.mo_ave_month, (SOTS_float_data_TEMP.flux_LD_corr_mo_ave/1000)*365,'-ok','LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LD_corr_mo_ave/1000)*365,'-o','color',cols(3,:),'LineWidth',2)

plot(SOTS_float_data_TEMP.mo_ave_month, (SOTS_float_data_TEMP.flux_LS_corr_mo_ave/1000)*365,'-ob','LineWidth',2)
plot(SOTS_float_data.mo_ave_month, (SOTS_float_data.flux_LS_corr_mo_ave/1000)*365,'-o','color',cols(5,:),'LineWidth',2)

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
plot(SOTS_float_data.time,SOTS_float_data.PSAL_20,'*b','LineWidth',1.8)
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_PSAL, '-.b')
ylabel('float salinity top 20m')

yyaxis right
plot(SOTS_float_data.time,SOTS_float_data.TEMP_20,'*r','LineWidth',1.8)
hold on
plot(mooring_data.xCO2_time, mooring_data.xCO2_SST, '-.r')
% legend('float', 'SOFS')
ylabel('float sw temperature top 20m')
xlim([datetime('10-12-2020', 'InputFormat','dd-MM-yyyy'), datetime('27-01-2022', 'InputFormat','dd-MM-yyyy')])
