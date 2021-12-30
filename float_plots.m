% comparison of float pH and CTD cast pH (based on DIC & Alk calculations)

load('SOTS_float_data.mat')

load('CTD_data.mat')

figure()
plot(SOTS_float_data.pH_LIR_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '^r', 'MarkerSize', 3)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), 'og', 'MarkerSize', 3)
plot(SOTS_float_data.pH_Williams_Deep(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), SOTS_float_data.pres(:,SOTS_float_data.time < datetime('18-12-2020','InputFormat','dd-MM-yyyy')), '+k', 'MarkerSize', 3)
plot(CTD_data.raw_data.pH(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
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
plot(CTD_data.raw_data.Alkalinity(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')), CTD_data.raw_data.Depth(CTD_data.raw_data.date < datetime('18-12-2020','InputFormat','dd-MM-yyyy') & CTD_data.raw_data.date > datetime('12-12-2020','InputFormat','dd-MM-yyyy')),'.b', 'MarkerSize', 12)
hold off
xlabel('Alk umol/kg')
ylabel('Depth dbar')
legend('Triangle float LIAR Alk','o float ES Alk','. CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
% ylim([0 100])


%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
subplot(3,2,1)
title('Red Triangle float pH LD, green = LD, cyan = WD, blue . CTD pH')
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(:,month(CTD_data.raw_data.date)==1), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==1),'.b','MarkerSize',10)
xlabel('pH total scale- Jan')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,2)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==2), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==2),'.b','MarkerSize',10)
xlabel('pH total scale - Feb')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,3)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==3), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==3),'.b','MarkerSize',10)
xlabel('pH total scale - Mar')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,4)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==4), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==4),'.b','MarkerSize',10)
xlabel('pH total scale - Apr')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,5)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==5), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==5),'.b','MarkerSize',10)
xlabel('pH total scale - May')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,6)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==6), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==6),'.b','MarkerSize',10)
xlabel('pH total scale - Jun')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off


figure()
subplot(3,2,1)
title('Red Triangle float pH LD, green = LD, cyan = WD, blue . CTD pH')
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==7), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==7),'.b','MarkerSize',10)
xlabel('pH total scale- July')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,2)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==8), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==8),'.b','MarkerSize',10)
xlabel('pH total scale - Aug')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,3)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==9), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==9),'.b','MarkerSize',10)
xlabel('pH total scale - Sept')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,4)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==10), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==10),'.b','MarkerSize',10)
xlabel('pH total scale - Oct')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,5)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==11), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==11),'.b','MarkerSize',10)
xlabel('pH total scale - Nov')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,6)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==12), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==12),'.b','MarkerSize',10)
xlabel('pH total scale - Dec')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 2000])
set(gca, 'YDir','reverse')
hold off

%%%%%%%%%
figure()
subplot(3,2,1)
title('Red Triangle float pH LD, green = LD, cyan = WD, blue . CTD pH')
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==1),SOTS_float_data.pres(:,month(SOTS_float_data.time)==1),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(:,month(CTD_data.raw_data.date)==1), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==1),'.b','MarkerSize',10)
xlabel('pH total scale- Jan')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,2)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==2),SOTS_float_data.pres(:,month(SOTS_float_data.time)==2),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==2), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==2),'.b','MarkerSize',10)
xlabel('pH total scale - Feb')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,3)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==3),SOTS_float_data.pres(:,month(SOTS_float_data.time)==3),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==3), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==3),'.b','MarkerSize',10)
xlabel('pH total scale - Mar')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,4)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==4),SOTS_float_data.pres(:,month(SOTS_float_data.time)==4),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==4), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==4),'.b','MarkerSize',10)
xlabel('pH total scale - Apr')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,5)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==5),SOTS_float_data.pres(:,month(SOTS_float_data.time)==5),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==5), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==5),'.b','MarkerSize',10)
xlabel('pH total scale - May')
ylabel('Depth dbar')
xlim([7.5 8.5])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,6)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==6),SOTS_float_data.pres(:,month(SOTS_float_data.time)==6),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==6), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==6),'.b','MarkerSize',10)
xlabel('pH total scale - Jun')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off


figure()
subplot(3,2,1)
title('Red Triangle float pH LD, green = LD, cyan = WD, blue . CTD pH')
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==7),SOTS_float_data.pres(:,month(SOTS_float_data.time)==7),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==7), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==7),'.b','MarkerSize',10)
xlabel('pH total scale- July')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,2)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==8),SOTS_float_data.pres(:,month(SOTS_float_data.time)==8),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==8), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==8),'.b','MarkerSize',10)
xlabel('pH total scale - Aug')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,3)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==9),SOTS_float_data.pres(:,month(SOTS_float_data.time)==9),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==9), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==9),'.b','MarkerSize',10)
xlabel('pH total scale - Sept')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,4)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==10),SOTS_float_data.pres(:,month(SOTS_float_data.time)==10),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==10), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==10),'.b','MarkerSize',10)
xlabel('pH total scale - Oct')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,5)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==11),SOTS_float_data.pres(:,month(SOTS_float_data.time)==11),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==11), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==11),'.b','MarkerSize',10)
xlabel('pH total scale - Nov')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off

subplot(3,2,6)
plot(SOTS_float_data.pH_LIR_Deep(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^r', 'MarkerSize',2)
hold on
plot(SOTS_float_data.pH_LIR_Shallow(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^g', 'MarkerSize',2)
plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==12),SOTS_float_data.pres(:,month(SOTS_float_data.time)==12),'^c', 'MarkerSize',2)
plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==12), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==12),'.b','MarkerSize',10)
xlabel('pH total scale - Dec')
ylabel('Depth dbar')
xlim([7.9 8.2])
ylim([0 100])
set(gca, 'YDir','reverse')
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
