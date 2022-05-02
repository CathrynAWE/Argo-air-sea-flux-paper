% comparison of float pH and CTD cast pH (based on DIC & Alk calculations)
clear all
load('SOTS_float_data.mat')
load('CTD_data.mat')
load('CTD_data_55S.mat')

%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot just the pH profiles for QC
% figure()
% plot(SOTS_float_data.pH_LIR_Deep(:,11), SOTS_float_data.pres(:,11),'or')
% hold on
% plot(SOTS_float_data.pH_LIR_Deep(:,11), SOTS_float_data.pres(:,11),'or')
% plot(SOTS_float_data.pH_LIR_Deep(:,23), SOTS_float_data.pres(:,23),'or')
% plot(SOTS_float_data.pH(:,10), SOTS_float_data.pres(:,10),'*b')
% plot(SOTS_float_data.pH(:,11), SOTS_float_data.pres(:,11),'*b')
% plot(SOTS_float_data.pH(:,23), SOTS_float_data.pres(:,23),'*b')
% 
% set(gca, 'YDir','reverse')
% xlabel('Depth')
% ylabel('pH total scale')
% ylim([0 100])
% 
%11, 23, 10 for SOTS float

%%% plot just the pH profiles for QC
% load('S55_float_data.mat')
% 
% figure()
% plot(S55_float_data.pH_LIR_Deep(:,1:3), S55_float_data.pres(:,1:3),'or')
% hold on
% plot(S55_float_data.pH_LIR_Deep(:,1:3), S55_float_data.pres(:,1:3),'or')
% hold on
% plot(S55_float_data.pH(:,1:3), S55_float_data.pres(:,1:3),'*b')
% 
% set(gca, 'YDir','reverse')
% xlabel('Depth')
% ylabel('pH total scale')
% ylim([0 100])
% 
% 1:3

%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.Depth <=1600), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')& CTD_data.raw_data.Depth <=1600),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
% legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% float pH vs SOLACE CTD cast pH
%%% this time with TS diagram and lat / lon of casts and float location
figure()
subplot(1,3,1)
plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
% plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
xlabel('pH total scale')
ylabel('Depth dbar')
title('SOTS float pH and SOLACE CTD casts')
legend('Triangle float LD','o float LS','+ float WD','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])

subplot(1,3,2)
plot(SOTS_float_data.TEMP(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.psal(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')

subplot(1,3,3)
plot(SOTS_float_data.lat(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.lon(:,SOTS_float_data.time < datetime('15-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data.raw_data.Latitude(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Longitude(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Latitude')
ylabel('Longitude')
legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')

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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOTS float estimated alkalinity vs CTD alkalinity
d =SOTS_float_data.pres <=100;
d(:,4:175)=0;

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
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.Depth <=1600), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('14-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')& CTD_data.raw_data.Depth <=1600),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
% legend('Triangle float','. CTD','Orientation','horizontal','Location','bestoutside')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% S55 float estimated alkalinity vs CTD alkalinity

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


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% all float pH vs CTD casts at SOTS
figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^g', 'MarkerSize',2)
%     plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^c', 'MarkerSize',2)
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==i), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==i),'.b','MarkerSize',10)
    if i ==1
        xlabel('pH total scale - Jan')
    elseif i ==2
        xlabel('pH total scale - Feb')
    elseif i==3
        xlabel('pH total scale - Mar')
    elseif i==4
        xlabel('pH total scale - Apr')
    elseif i==5
        xlabel('pH total scale - May')
    elseif i==6
        xlabel('pH total scale - June')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
    
end

figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^g', 'MarkerSize',2)
%     plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^c', 'MarkerSize',2)
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==(i+6)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==(i+6)),'.b','MarkerSize',10)
    if i ==1
        xlabel('pH total scale - July')
    elseif i ==2
        xlabel('pH total scale - Aug')
    elseif i==3
        xlabel('pH total scale - Sep')
    elseif i==4
        xlabel('pH total scale - Oct')
    elseif i==5
        xlabel('pH total scale - Nov')
    elseif i==6
        xlabel('pH total scale - Dec')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% all float pH vs CTD casts at SOTS with temp / sal trace
%%% Mar, Apr, July, Aug, Oct, Nov, Dec
figure()
n=[3,3,4,4,7,7];
for i = 1:6
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LD_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^m', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LS_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^y', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',12)
        ylabel('Depth dbar')
        xlim([7.95 8.1])
        ylim([0 100])
        set(gca, 'YDir','reverse')
        if i ==1
            xlabel('pH total scale - Mar')
        elseif i ==3
            xlabel('pH total scale - Apr')
        elseif i==5
            xlabel('pH total scale - July')
        end
        title('SOTS float')
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
        title('SOTS float')
        hold off
        
    end
end

%%% Mar, Apr, July, Aug, Oct, Nov, Dec
figure()
n=[8,8,10,10,11,11];
for i = 1:6
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LD_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^m', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LS_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^y', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',12)
        ylabel('Depth dbar')
        xlim([7.95 8.1])
        ylim([0 100])
        set(gca, 'YDir','reverse')
        if i ==1
            xlabel('pH total scale - Aug')
        elseif i ==3
            xlabel('pH total scale - Oct')
        elseif i==5
            xlabel('pH total scale - Nov')
        end
        title('SOTS float')
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
        title('SOTS float')
        hold off
 
    end
end


figure()
n=[12,12];
for i = 1:2
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LD_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^m', 'MarkerSize',2)
        plot(SOTS_float_data.pH_LS_corr(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^y', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',12)
        ylabel('Depth dbar')
        xlim([7.95 8.1])
        ylim([0 100])
        set(gca, 'YDir','reverse')
%         if i ==1
            xlabel('pH total scale - Dec')
%         elseif i ==3
%             xlabel('pH total scale - Oct')
%         elseif i==5
%             xlabel('pH total scale - Nov')
%         end
        title('SOTS float')
%     
    elseif rem(i,2)==0
        subplot(3,2,i)
subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
        title('SOTS float')
        hold off

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% all float pH vs CTD casts at SOTS with lat / lon of both
%%% Mar, Apr, July, Aug, Oct, Nov, Dec
figure()
n=[3,3,4,4,7,7];
for i = 1:6
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Depth dbar')
        xlim([7.9 8.2])
        ylim([0 100])
        set(gca, 'YDir','reverse')
        if i ==1
            xlabel('pH total scale - Mar')
        elseif i ==3
            xlabel('pH total scale - Apr')
        elseif i==5
            xlabel('pH total scale - July')
        end
        title('SOTS float')
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.lat(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.lon(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(CTD_data.raw_data.Latitude(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Longitude(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Latitude')
        xlabel('Longitude')
        title('SOTS float')
        hold off
        
    end
end

%%% Mar, Apr, July, Aug, Oct, Nov, Dec
figure()
n=[8,8,10,10,11,11];
for i = 1:6
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Depth dbar')
        xlim([7.9 8.2])
        ylim([0 100])
        set(gca, 'YDir','reverse')
        if i ==1
            xlabel('pH total scale - Aug')
        elseif i ==3
            xlabel('pH total scale - Oct')
        elseif i==5
            xlabel('pH total scale - Nov')
        end
        title('SOTS float')
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.lat(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.lon(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(CTD_data.raw_data.Latitude(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Longitude(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Latitude')
        xlabel('Longitude')
        title('SOTS float')
        hold off
 
    end
end


figure()
n=[12,12];
for i = 1:2
    if rem(i,2)~=0
        subplot(3,2,i)
        plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^g', 'MarkerSize',2)
        plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==n(i)),'^c', 'MarkerSize',2)
        plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==n(i)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Depth dbar')
        xlim([7.9 8.2])
        ylim([0 100])
        set(gca, 'YDir','reverse')
%         if i ==1
            xlabel('pH total scale - Dec')
%         elseif i ==3
%             xlabel('pH total scale - Oct')
%         elseif i==5
%             xlabel('pH total scale - Nov')
%         end
       title('SOTS float')
    elseif rem(i,2)==0
        subplot(3,2,i)
        subplot(3,2,i)
        plot(SOTS_float_data.lat(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.lon(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot(CTD_data.raw_data.Latitude(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Longitude(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('Latitude')
        xlabel('Longitude')
        title('SOTS float')
        hold off

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% all float alkalinity vs CTD casts at SOTS
figure()
for i = 1:6
    subplot(3,2,i)
        plot(SOTS_float_data.Alk_LIAR(:,month(SOTS_float_data.time)==i),SOTS_float_data.Alk_pres(:,month(SOTS_float_data.time)==i),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.Alk_ES(:,month(SOTS_float_data.time)==i),SOTS_float_data.Alk_pres(:,month(SOTS_float_data.time)==i),'^g', 'MarkerSize',2)
    plot(CTD_data.raw_data.Alkalinity(CTD_data.raw_data.month==i), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==i),'.b','MarkerSize',10)
    if i ==1
        xlabel('Alkalinity - Jan')
    elseif i ==2
        xlabel('Alkalinity - Feb')
    elseif i==3
        xlabel('Alkalinity - Mar')
    elseif i==4
        xlabel('Alkalinity - Apr')
    elseif i==5
        xlabel('Alkalinity - May')
    elseif i==6
        xlabel('Alkalinity - June')
    end
    title('SOTS float - Red Triangle float Alk LIAR, green = Alk ES, blue . CTD pH')
    ylabel('Depth dbar')
    % xlim([7.9 8.2])
    ylim([0 100])
    set(gca, 'YDir','reverse')
    
end

figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data.Alk_LIAR(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.Alk_pres(:,month(SOTS_float_data.time)==(i+6)),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.Alk_ES(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.Alk_pres(:,month(SOTS_float_data.time)==(i+6)),'^g', 'MarkerSize',2)
    plot(CTD_data.raw_data.Alkalinity(CTD_data.raw_data.month==(i+6)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==(i+6)),'.b','MarkerSize',10)
    if i ==1
        xlabel('Alkalinity - July')
    elseif i ==2
        xlabel('Alkalinity - Aug')
    elseif i==3
        xlabel('Alkalinity - Sep')
    elseif i==4
        xlabel('Alkalinity - Oct')
    elseif i==5
        xlabel('Alkalinity - Nov')
    elseif i==6
        xlabel('Alkalinity - Dec')
    end
    ylabel('Depth dbar')
    title('SOTS float - Red Triangle float Alk LIAR, green = Alk ES, blue . CTD pH')
    % xlim([7.9 8.2])
    ylim([0 100])
    set(gca, 'YDir','reverse')
    
end


%%%%%%%%%%%%%%%%
% monthly averages of the top 20m

figure()
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_LD_20_mo_ave,'^r', 'MarkerSize',4)
hold on
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_LS_20_mo_ave,'^g', 'MarkerSize',4)
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_WD_20_mo_ave,'^k', 'MarkerSize',4)
plot(CTD_data.pH_monthly_ave_20.month,CTD_data.pH_monthly_ave_20.pH(CTD_data.pH_monthly_ave_20.pH~=0),'.b','MarkerSize',10)
legend('float pH - LD','float pH - LS','float pH - WD','CTD','Orientation',...
    'horizontal','Location','bestoutside','FontSize',12)
xlabel('Month')
ylabel('pH total scale')
title('SOTS float monthly ave top 20m only')
hold off


figure()
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.Alk_LIAR_20_mo_ave,'^r', 'MarkerSize',4)
hold on
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.Alk_ES_20_mo_ave,'^g', 'MarkerSize',4)
plot(CTD_data.Alk_monthly_ave_20.month,CTD_data.Alk_monthly_ave_20.Alk(CTD_data.Alk_monthly_ave_20.Alk~=0),'.b','MarkerSize',10)
legend('float Alk - LIAR','float Alk - ES','CTD','Orientation',...
    'horizontal','Location','bestoutside','FontSize',12)
xlabel('Month')
xlim([0 13])
ylabel('Alkalinity umol/kg')
title('SOTS float monthly ave top 20m only')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% all float pCO2 data as calculated from the various pH corrections and
%%%%% Alkalinity estimates
load('mooring_data.mat')
load('SOTS_float_data.mat')

figure()
plot(SOTS_float_data.time,SOTS_float_data.pCO2_uatm,'^k','MarkerSize',4)
hold on
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D,'.r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S,'.g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D,'.c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_corr,'+r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_corr,'+g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_corr,'+c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_ES,'*r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_ES,'*g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_ES,'*c','MarkerSize',6)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% same for S55 float

%%%%% all float pCO2 data as calculated from the various pH corrections and
%%%%% Alkalinity estimates
clear all
load('S55_float_data.mat')

figure()
% plot(S55_float_data.time,S55_float_data.pCO2_uatm,'^k','MarkerSize',4)
hold on
% plot(S55_float_data.time,S55_float_data.pCO2_L_D,'.r','MarkerSize',6)
% plot(S55_float_data.time,S55_float_data.pCO2_L_S,'.g','MarkerSize',6)
% plot(S55_float_data.time,S55_float_data.pCO2_W_D,'.c','MarkerSize',6)
plot(S55_float_data.time,S55_float_data.pCO2_L_D_corr,'+r','MarkerSize',6)
plot(S55_float_data.time,S55_float_data.pCO2_L_S_corr,'+g','MarkerSize',6)
plot(S55_float_data.time,S55_float_data.pCO2_W_D_corr,'+c','MarkerSize',6)
% plot(S55_float_data.time,S55_float_data.pCO2_L_D_ES,'*r','MarkerSize',6)
% plot(S55_float_data.time,S55_float_data.pCO2_L_S_ES,'*g','MarkerSize',6)
% plot(S55_float_data.time,S55_float_data.pCO2_W_D_ES,'*c','MarkerSize',6)
hold off
xlabel('Time')
xlim([datetime('01-12-2020','inputFormat','dd-MM-yyyy') datetime('15-01-2022','inputFormat','dd-MM-yyyy')])
ylabel('pCO_2 uatm')
legend('Cape Grim atm pCO_2','float LD-LIAR','float LS-LIAR',...
    'float WD-LIAR','float LD corr','float LS corr','float WD corr',...
    'Orientation','vertical','Location','westoutside','FontSize',8)
title('S55 float')
% 'float pH-LD-ES', 'float pH-LS-ES','float pH-WD-ES',




%%%%%%%%%%%%%%%%%%%
%%%% mooring wsp vs ERA5 and NCEP wsp
clear all
load('SOTS_float_data.mat')
load('mooring_data.mat')
load('NCEP_SOTS_float.mat')

figure()
plot(SOTS_float_data.time, SOTS_float_data.wsp,'^r','MarkerSize',4)
hold on
plot(mooring_data.wsp_pCO2_time,mooring_data.wsp_pCO2,'*b','MarkerSize',2)
plot(SOTS_float_data.time, float_ncep.wnd_ncep_f,'ok')
hold off
xlabel('Time')
ylabel('Windspeed m^-2')
title('SOTS float')
legend('ERA5 wsp','mooring wsp','NCEP wsp','Orientation','horizontal','Location','southoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% mooring wsp vs ERA5 and NCEP wsp
clear all

load('mooring_data.mat')
load('S55_float_data.mat')
load('NCEP_S55_float.mat')

figure()
plot(S55_float_data.time, S55_float_data.wsp,'^r','MarkerSize',4)
hold on
plot(S55_float_data.time, float_ncep.wnd_ncep_f,'ok')
hold off
xlabel('Time')
ylabel('Windspeed m^-2')
title('55S float')
legend('ERA5 wsp','NCEP wsp','Orientation','horizontal','Location','southoutside')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% mooring monthly average fluxes vs float monthly average fluxes
% figure()
% plot([1:12],mooring_data.pCO2_2021_monthly_flux,'-r')
% hold on
% plot([1:12],mooring_data.pCO2_2020_monthly_flux,'-b')
% plot(SOTS_float_data.mo_ave_month, SOTS_float_data.flux_LD_mo_ave,'--g')
% plot(SOTS_float_data.mo_ave_month, SOTS_float_data.flux_LS_mo_ave,'--c')
% plot(SOTS_float_data.mo_ave_month, SOTS_float_data.flux_WD_mo_ave,'--k')
% hold off
% legend('mooring 2021', 'mooring 2020', 'float LD 2020', 'float LS 2020',...
%     'float WD 2020')
% xlabel('Month')
% ylabel('air sea flux mol m^-2 d^-^1')
% xticks([0:13])
% xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
%     'Sep', 'Oct', 'Nov', 'Dec',''}) 
% xlim([0 13])


%%%%%%%%%%%%%%%%
%%%% ERA5 windspeed, air temp and mean sea level pressure vs NCEP
%%%% reanalysis data
clear all
load('SOTS_float_data.mat')
load('NCEP_SOTS_float.mat')
figure()
yyaxis left
plot(SOTS_float_data.time, SOTS_float_data.wsp,'og')
hold on
plot(SOTS_float_data.time, float_ncep.wnd_ncep_f,'ok')
% plot(mooring_data.wsp_time,mooring_data.wsp,'ob')
plot(SOTS_float_data.time, SOTS_float_data.t_2,'+r')
plot(SOTS_float_data.time, float_ncep.airT_ncep_f,'+b')

hold off
ylabel('windspeed m s^-^1 air temp C')
title('SOTS float')

yyaxis right
plot(SOTS_float_data.time, SOTS_float_data.msp, '*r')
hold on
plot(SOTS_float_data.time, float_ncep.mslp_ncep_f,'*b')
hold off
ylabel('mean sea level pressure Pa')
xlabel('Time')
legend('wsp ERA5', 'wsp NCEP', 'air T ERA5', 'air T NCEP', 'mslp ERA5',...
    'mslp NCEP')
xlim([datetime('01-01-2021','InputFormat','dd-MM-yyyy') datetime('31-12-2021','InputFormat','dd-MM-yyyy')]) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% same for 55S float
%%%% ERA5 windspeed, air temp and mean sea level pressure vs NCEP
%%%% reanalysis data
clear all
load('S55_float_data.mat')
load('NCEP_S55_float.mat')

figure()
yyaxis left
plot(S55_float_data.time, S55_float_data.wsp,'og')
hold on
plot(S55_float_data.time, float_ncep.wnd_ncep_f,'ok')
plot(S55_float_data.time, S55_float_data.t_2,'+r')
plot(S55_float_data.time, float_ncep.airT_ncep_f,'+b')
hold off
ylabel('windspeed m s^-^1 air temp C')
title('55S float')

yyaxis right
plot(S55_float_data.time, S55_float_data.msp, '*r')
hold on
plot(S55_float_data.time, float_ncep.mslp_ncep_f,'*b')
hold off
ylabel('mean sea level pressure Pa')
xlabel('Time')
legend('wsp ERA5', 'wsp NCEP', 'air T ERA5', 'air T NCEP', 'mslp ERA5',...
    'mslp NCEP')
% xlim([datetime('01-01-2021','InputFormat','dd-MM-yyyy') datetime('31-12-2021','InputFormat','dd-MM-yyyy')]) 