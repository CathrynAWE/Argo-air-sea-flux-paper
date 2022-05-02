float_ID = '5906623' % SOTS float
% float_ID = '5906624' % 55S float
%%%%% change structure name!

addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0\seawater_ver3_0'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper'

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)


% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name];

% read necessary variables from float profile
lat = ncread(file,'LATITUDE');
lon = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
JULD = ncread(file,'JULD');
pres = ncread(file,'PRES'); %dbar
temp = ncread(file, 'TEMP_ADJUSTED');% Celsius
temp_QC = ncread(file, 'TEMP_ADJUSTED_QC');
oxy = ncread(file,'DOXY_ADJUSTED');% umol/kg
oxy_QC = ncread(file,'DOXY_ADJUSTED_QC');
NO3 = ncread(file,'NITRATE');
pH_adjusted = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH = ncread(file,'PH_IN_SITU_TOTAL');
% pH_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
N = ncread(file, 'NITRATE_ADJUSTED');
S = ncread(file,'PSAL');
S_QC = ncread(file,'PSAL_QC');
Cycle = ncread(file,'CYCLE_NUMBER');
fs = size(pH);


for i=1:fs(1,2)
    idx_950(:,i) = pres(:,i)>900 & pres(:,i)<1000;
    idx_1500(:,i) = pres(:,i)>1450 & pres(:,i)<1550;
end

figure()
hold on
for i =1:fs(1,2)
    if sum(idx_950(:,i))>0
        yyaxis left
        plot(temp(idx_950(:,i),i),S(idx_950(:,i),i),'ok')
        ylabel('salinity')
        yyaxis right
        plot(temp(idx_950(:,i),i),oxy(idx_950(:,i),i),'*k')
        ylabel('dissolved oxygen - umol/kg')
        xlabel('temperature')
    end
    if sum(idx_1500(:,i))>0
        yyaxis left
        plot(temp(idx_1500(:,i),i),S(idx_1500(:,i),i),'or')
        yyaxis right
        plot(temp(idx_1500(:,i),i),oxy(idx_1500(:,i),i),'*r')
     end
    
end

