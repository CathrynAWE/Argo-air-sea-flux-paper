% comparison of float pH and CTD cast pH (based on DIC & Alk calculations)

load('SOTS_float_data.mat')

load('CTD_data.mat')

figure()
subplot(2,1,1)
plot(SOTS_float_data.time,SOTS_float_data.pH_adj_minus,'or', 'MarkerSize',2)
hold on
plot(CTD_data.date, CTD_data.pH,'-b')
xlabel('Time')
ylabel('pH in situ total')
legend('float', 'CTD','Orientation', 'horizontal','Location','bestoutside')
hold off

subplot(2,1,2)
yyaxis left
plot(SOTS_float_data.time, SOTS_float_data.lat, 'o','MarkerSize',2)
hold on
plot(CTD_data.date, CTD_data.Latitude, '-')
ylabel('Latitude')
hold off

yyaxis right
plot(SOTS_float_data.time, SOTS_float_data.lon, 'o','MarkerSize',2)
hold on
plot(CTD_data.date, CTD_data.Longitude, '-')
ylabel('Longitude')
hold off

figure()
plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data.pH(CTD_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.Depth(CTD_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
legend('Triangle float LD','o float LS','+ float WD','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])


figure()
plot(SOTS_float_data.Alk_LIAR(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.Alk_ES(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.Alk_pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(CTD_data.Alkalinity(CTD_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.Depth(CTD_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Alk umol/kg')
ylabel('Depth dbar')
legend('Triangle float LIAR Alk','o float ES Alk','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 100])