figure()
subplot(2,1,1)
title('SOTS 2020 & 2021 UW data')
yyaxis left
plot(time_SOTS2021,lat_SOTS2021,'-r')
hold on
plot(time_SOTS2020,lat_SOTS2020,'--r')
hold off
ylabel('Latitude')
yyaxis right
plot(time_SOTS2021,lon_SOTS2021,'-b')
hold on
plot(time_SOTS2020,lon_SOTS2020,'--b')
hold off
ylabel('Longitude')
xlabel('Time')

subplot(2,1,2)
yyaxis left
plot(time_SOTS2021,F_CO2_SOTS2021,'ok','MarkerSize',2)
hold on
plot(time_SOTS2020,F_CO2_SOTS2020,'.b','MarkerSize',4)
hold off
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
yyaxis right
plot(time_SOTS2021,T_SOTS2021,'-b')
hold on
plot(time_SOTS2020,T_SOTS2020,'--b')
hold off
ylabel('SST')
xlabel('Time')