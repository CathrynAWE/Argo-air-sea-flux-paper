% comparison of float pH and CTD cast pH (based on DIC & Alk calculations)
clear all
load('S55_comp_float_data.mat')


%%%%%%%%
%%% plot the various corrected pH against one another
figure()
plot(S55_comp_float_data.pH_adj,S55_comp_float_data.pres,'or')
hold on
plot(S55_comp_float_data.pH_LIR_Deep, S55_comp_float_data.pres, '*b')
plot(S55_comp_float_data.pH_LIR_Shallow, S55_comp_float_data.pres,'+g')
% plot(S55_comp_float_data.pH_LD_corr, S55_comp_float_data.pres,'.b')
% plot(S55_comp_float_data.pH_LS_corr,S55_comp_float_data.pres,'^g')
hold off
set(gca,'YDir','reverse')
xlabel('pressure')
ylabel('pH total scale')
legend('red pH adj', 'blue pH LD', 'green pH LS',  'Orientation',...
    'horizontal','Location','southoutside')
%'blue pH LD corr','green pH LS corr',
xlim([7.8 8.05])


%%%%%%
%%% plot the different pCO2_sw from the different pH against one another
figure()
plot(S55_comp_float_data.time, S55_comp_float_data.pCO2_pH_adj,'or')
hold on
plot(S55_comp_float_data.time, S55_comp_float_data.pCO2_L_D,'*b')
plot(S55_comp_float_data.time, S55_comp_float_data.pCO2_L_S,'+g')
plot(S55_comp_float_data.time, S55_comp_float_data.pCO2_L_D_corr,'.b')
plot(S55_comp_float_data.time, S55_comp_float_data.pCO2_L_S_corr,'.g')
hold off
xlabel('Time')
ylabel('pCO_2 sw')
legend('red pCO2 pH adj', 'blue pCO2 pH LD', 'green pCO2 pH LS',...
     'blue pCO2 pH LD corr', 'green pCO2 pH LS corr','Orientation',...
    'horizontal','Location','southoutside')
ylim([350 450])


%%%%%
%%% plot the different air sea fluxes vs S55 float from 2021
load('S55_float_data.mat')

figure()
plot(S55_comp_float_data.time, S55_comp_float_data.flux_pH_adj,'or')
hold on
plot(S55_comp_float_data.time, S55_comp_float_data.flux_L_D,'*b')
plot(S55_comp_float_data.time, S55_comp_float_data.flux_L_S,'+g')
plot(S55_float_data.time, S55_float_data.flux_L_D,'.b')
plot(S55_float_data.time, S55_float_data.flux_L_S,'.g')
hold off
xlabel('Time')
ylabel('air sea flux mmol m^-^2 d^-^1')
legend('red flux pH adj', 'blue flux pH LD', 'green flux pH LS',...
    'blue . flux pH LD S55', 'green . flux pH LS S55','','Orientation',...
    'horizontal','Location','southoutside')

%%%%%
%%% plot the different air sea fluxes vs S55 float as monthly means
figure()
plot([1:12],S55_comp_float_data.y2016.flux_pH_adj_mo_ave,'*r')
hold on
plot([1:12],S55_comp_float_data.y2016.flux_LD_mo_ave,'*b')
plot([1:12],S55_comp_float_data.y2016.flux_LS_mo_ave,'*g')
plot([1:12],S55_comp_float_data.y2017.flux_pH_adj_mo_ave,'or')
plot([1:12],S55_comp_float_data.y2017.flux_LD_mo_ave,'ob')
plot([1:12],S55_comp_float_data.y2017.flux_LS_mo_ave,'og')
plot(S55_comp_float_data.y2018.mo_ave_month,S55_comp_float_data.y2018.flux_pH_adj_mo_ave,'.r')
plot(S55_comp_float_data.y2018.mo_ave_month,S55_comp_float_data.y2018.flux_LD_mo_ave,'.b')
plot(S55_comp_float_data.y2018.mo_ave_month,S55_comp_float_data.y2018.flux_LS_mo_ave,'.g')
plot([1:12],S55_float_data.flux_LD_mo_ave,'+c')
plot([1:12],S55_float_data.flux_LS_mo_ave,'+m')
hold off
ylabel('air sea flux mmol m^-^2 d^-^1')
legend('red * pH adj 2016', 'blue * pH LD 2016', 'green * pH LS 2016',...
'red o pH adj 2017', 'blue o pH LD 2017', 'green o pH LS 2017',...
'red . pH adj 2018', 'blue . pH LD 2018', 'green . pH LS 2018',...
'cyan + pH LD S55', 'magenta + pH LS S55','','Orientation',...
    'vertical','Location','westoutside')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)
