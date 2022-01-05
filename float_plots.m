% comparison of float pH and CTD cast pH (based on DIC & Alk calculations)

load('SOTS_float_data.mat')

load('CTD_data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%
%%% float pH vs SOLACE CTD cast pH

figure()
title('float pH and SOLACE CTD casts')
plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
legend('Triangle float LD','o float LS','+ float WD','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])

%%% this time with TS diagram
figure()
subplot(1,2,1)
title('float pH and SOLACE CTD casts')
plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
legend('Triangle float LD','o float LS','+ float WD','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])

subplot(1,2,2)
plot(SOTS_float_data.TEMP(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.psal(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(CTD_data.raw_data.T_insitu(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Salinity_CTD(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('09-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Temp')
ylabel('PSAL')
legend('Triangle float temp','. CTD temp','Orientation','horizontal','Location','bestoutside')



%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
plot(SOTS_float_data.Alk_LIAR(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.Alk_ES(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(CTD_data.raw_data.Alkalinity(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Alk umol/kg')
ylabel('Depth dbar')
legend('Triangle float LIAR Alk','o float ES Alk','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 100])


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% all float pH vs CTD casts at SOTS
figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^g', 'MarkerSize',2)
    plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==i),SOTS_float_data.pres(:,month(SOTS_float_data.time)==i),'^c', 'MarkerSize',2)
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
    xlim([7.9 8.2])
    ylim([0 100])
    set(gca, 'YDir','reverse')
    
end

figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^g', 'MarkerSize',2)
    plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^c', 'MarkerSize',2)
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
    xlim([7.9 8.2])
    ylim([0 100])
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
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
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
    
    elseif rem(i,2)==0
        subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
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
%     
    elseif rem(i,2)==0
        subplot(3,2,i)
subplot(3,2,i)
        plot(SOTS_float_data.TEMP(:,month(SOTS_float_data.time)==n(i)),SOTS_float_data.psal(:,month(SOTS_float_data.time)==n(i)),'^r', 'MarkerSize',2)
        hold on
        plot( CTD_data.raw_data.T_insitu(month(CTD_data.raw_data.date)==n(i)),CTD_data.raw_data.Salinity_CTD(month(CTD_data.raw_data.date)==n(i)),'.b','MarkerSize',10)
        ylabel('PSAL')
        xlabel('TEMP')
        hold off

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% all float alkalinity vs CTD casts at SOTS
figure()
for i = 1:6
    subplot(3,2,i)
    title('Red Triangle float Alk LIAR, green = Alk ES, blue . CTD pH')
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
    ylabel('Depth dbar')
    % xlim([7.9 8.2])
    ylim([0 100])
    set(gca, 'YDir','reverse')
    
end

figure()
for i = 1:6
    subplot(3,2,i)
    title('Red Triangle float Alk LIAR, green = Alk ES, blue . CTD pH')
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
    % xlim([7.9 8.2])
    ylim([0 100])
    set(gca, 'YDir','reverse')
    
end


%%%%%%%%%%%%%%%%
% monthly averages of the top 20m

figure()
title('monthly ave top 20m only')
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_LD_20_mo_ave,'^r', 'MarkerSize',4)
hold on
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_LS_20_mo_ave,'^g', 'MarkerSize',4)
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.pH_WD_20_mo_ave,'^k', 'MarkerSize',4)
plot(CTD_data.pH_monthly_ave_20.month,CTD_data.pH_monthly_ave_20.pH(CTD_data.pH_monthly_ave_20.pH~=0),'.b','MarkerSize',10)
legend('float pH - LD','float pH - LS','float pH - WD','CTD','Orientation',...
    'horizontal','Location','bestoutside','FontSize',12)
xlabel('Month')
ylabel('pH total scale')
hold off


figure()
title('monthly ave top 20m only')
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.Alk_LIAR_20_mo_ave,'^r', 'MarkerSize',4)
hold on
plot(SOTS_float_data.mo_ave_month,SOTS_float_data.Alk_ES_20_mo_ave,'^g', 'MarkerSize',4)
plot(CTD_data.Alk_monthly_ave_20.month,CTD_data.Alk_monthly_ave_20.Alk(CTD_data.Alk_monthly_ave_20.Alk~=0),'.b','MarkerSize',10)
legend('float Alk - LIAR','float Alk - ES','CTD','Orientation',...
    'horizontal','Location','bestoutside','FontSize',12)
xlabel('Month')
xlim([0 13])
ylabel('Alkalinity umol/kg')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% all float pCO2 data as calculated from the various pH corrections and
%%%%% Alkalinity estimates

figure()
plot(SOTS_float_data.time,SOTS_float_data.pCO2_atm,'^k','MarkerSize',4)
hold on
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D,'.r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S,'.g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D,'.c','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_D_ES,'*r','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_L_S_ES,'*g','MarkerSize',6)
plot(SOTS_float_data.time,SOTS_float_data.pCO2_W_D_ES,'*c','MarkerSize',6)
hold off
xlabel('Time')
xlim([datetime('01-12-2020','inputFormat','dd-MM-yyyy') datetime('15-01-2022','inputFormat','dd-MM-yyyy')])
ylabel('pCO_2 atm')
legend('Cape Grim atm pCO_2','float pH - LD - LIAR','float pH - LS - LIAR',...
    'float pH - WD - LIAR', 'float pH - LD - ES', 'float pH - LS - ES',...
    'float pH - WD - ES','Orientation','horizontal','Location','southoutside',...
    'FontSize',8)


%%%%%%%%%%%%%%%%%%%
%%%% mooring wsp vs ERA5 wsp
load('mooring_data.mat')

figure()
plot(SOTS_float_data.time, SOTS_float_data.wsp,'^r','MarkerSize',4)
hold on
plot(mooring_data.wsp_pCO2_time,mooring_data.wsp_pCO2,'*b','MarkerSize',2)
hold off
xlabel('Time')
ylabel('Windspeed m^-2')
legend('ERA5 wsp','mooring wsp','Orientation','horizontal','Location','southoutside')
