% this script calculates CTD cast pH to compare with 55S float pH
clear all

CTD_data_55S=[];

file = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\CTD_data\IN2010_2020_carbon results.xlsx';

data = readtable(file, 'Sheet','55S','ReadVariableNames',true);

CTD_data_55S.raw_data = data;

CTD_data_55S.raw_data.date = datetime(data.Date_Time, 'InputFormat','dd/MM/yyyy HH:mm');

CTD_data_55S.raw_data.month = month(CTD_data_55S.raw_data.date);


%Si ave from RAS SOFS8 2.8 umol/kg
%PO4 ave from RAS SOFS8 0.9
%NH4 set to 2 umol/kg
%H2S set to 0
% CTD salinity
[DATA, HEADERS, NICEHEADERS]  = CO2SYS(CTD_data_55S.raw_data.Alkalinity,CTD_data_55S.raw_data.TCO2,1,2,CTD_data_55S.raw_data.Salinity_CTD,CTD_data_55S.raw_data.T_insitu,CTD_data_55S.raw_data.T_insitu,CTD_data_55S.raw_data.Depth,CTD_data_55S.raw_data.Depth,2.8,0.9,2,0,1,10,1,2,2);
CTD_data_55S.raw_data.pH = DATA(:,21); %total scale
% change the -999 of pH column to NaN
CTD_data_55S.raw_data.pH(CTD_data_55S.raw_data.pH==-999)=NaN;


% sort it so that I can use accumarray and get the order of the months
CTD_data_55S.raw_data_sorted = sortrows(CTD_data_55S.raw_data,16);


% just the top 20m
CTD_data_55S.raw_data_sorted_20=CTD_data_55S.raw_data_sorted(CTD_data_55S.raw_data_sorted.Depth<=20,:);


CTD_data_55S.Alk_monthly_ave_20.Alk = accumarray(CTD_data_55S.raw_data_sorted_20.month, CTD_data_55S.raw_data_sorted_20.Alkalinity,[],@(x)mean(x,'omitnan'));
CTD_data_55S.Alk_monthly_ave_20.month = unique(CTD_data_55S.raw_data_sorted_20.month);

CTD_data_55S.DIC_monthly_ave_20.DIC = accumarray(CTD_data_55S.raw_data_sorted_20.month, CTD_data_55S.raw_data_sorted_20.TCO2,[],@(x)mean(x,'omitnan'));
CTD_data_55S.DIC_monthly_ave_20.month = unique(CTD_data_55S.raw_data_sorted_20.month);

CTD_data_55S.pH_monthly_ave_20.pH = accumarray(CTD_data_55S.raw_data_sorted_20.month, CTD_data_55S.raw_data_sorted_20.pH,[],@(x)mean(x,'omitnan'));
CTD_data_55S.pH_monthly_ave_20.month = unique(CTD_data_55S.raw_data_sorted_20.month);



clearvars -except CTD_data_55S

path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper';
cd(path)

save('CTD_data_55S.mat')

