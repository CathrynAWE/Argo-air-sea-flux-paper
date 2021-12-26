% this script calculates CTD cast pH to compare with float pH

file = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\CTD_data\IN2020_V08 V09_carbon results.xlsx';

CTD_data = readtable(file, 'Sheet','Results','ReadVariableNames',true);

CTD_data.date = datetime(CTD_data.Date_Time, 'InputFormat','dd/MM/yyyy HH:mm');


%Si ave from RAS SOFS8 2.8 umol/kg
%PO4 ave from RAS SOFS8 0.9
%NH4 set to 2 umol/kg
%H2S set to 0
% SOMMA salinity

[DATA, HEADERS, NICEHEADERS]  = CO2SYS(CTD_data.Alkalinity,CTD_data.TCO2,1,2,CTD_data.Salinity_CTD,CTD_data.T_insitu,CTD_data.T_insitu,CTD_data.Depth,CTD_data.Depth,2.8,0.9,2,0,1,10,1,2,2);
CTD_data.pH = DATA(:,21); %total scale
% CTD_data.pH_2 = DATA(:,43); %total scale, these two are the same column
% in DATA

clearvars -except CTD_data

path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper';
cd(path)

save('CTD_data.mat')

