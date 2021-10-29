% convert float salinity into alkalinity and then airsea flux
% now go find the data
% float IDs, choose one: 5906623 (SOTS), 5906624 (55S)
clear all

float_ID = '5906624'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\seawater_ver3_0\seawater_ver3_0'
search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name]


lat = ncread(file,'LATITUDE');
lon = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
pres = ncread(file,'PRES');
temp = ncread(file, 'TEMP_ADJUSTED');
oxy = ncread(file,'DOXY_ADJUSTED');
NO3 = ncread(file,'NITRATE');
pH = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
N = ncread(file, 'NITRATE_ADJUSTED');
T = ncread(file,'TEMP');
S = ncread(file,'PSAL');
fs = size(pH);

Alk_LIAR=[];
Alk_LIAR_20=[];
pH_20 =[];
S_20 = [];
temp_20 = [];
pres_20 =[];
for i = 1:fs(2)
%for i =1
    % mask out the NaNs in the float data, for the Alklinity script
    msk = ~isnan(oxy(:,i));
    lon_c = zeros(size(oxy(msk,i),1),1);
    lon_c(:) = lon(i,1);
    lat_c = zeros(size(oxy(msk,i),1),1);
    lat_c(:) = lat(i,1);
    Coor = [lon_c lat_c pres(msk,i)];

    Meas = [S(msk,i) oxy(msk,i) T(msk,i)];
    
    % define the measurement IDs (salinity, oxygen, temperature)
    MeasID = [1 6 7];

    % run the LIAR code to create total alkalinity     
    Alk_LIAR(1:size(oxy(msk,i),1),end+1) = LIAR(Coor, Meas, MeasID);
    
    
    % now we just want the top 20m average for the CO2SYS calcluations
    idx = pres(msk,i)<20;
    % a complicated way of getting an index vector for the accumarray
    % function
    c=ones(length(idx),1);
    c=c(idx);
    b=zeros(length(Alk_LIAR(:,1)),1);
    b(c~=0)=c;
    b=b+1;
    
    e=zeros(length(pres(msk,i)),1);
    e(c~=0)=c;
    e=e+1;
    
    % we only want the first value of accumarray (i.e. the first 20m)  
    a = accumarray(b,Alk_LIAR(:,i),[],@(x)mean(x,'omitnan'));
    Alk_LIAR_20(end+1,1) = a(1);
    p = accumarray(e,pH(msk,i),[],@(x)mean(x,'omitnan'));
    pH_20(end+1,1) = p(1);
    s = accumarray(e,S(msk,i),[],@(x)mean(x,'omitnan'));
    S_20(end+1,1) = s(1);
    pp = accumarray(e,pres(msk,i),[],@(x)mean(x,'omitnan'));
    pres_20(end+1,1) = pp(1);
    t = accumarray(e,temp(msk,i),[],@(x)mean(x,'omitnan'));
    temp_20(end+1,1) = t(1);
    
end
Alk_LIAR(Alk_LIAR==0)=NaN;

%Si ave from RAS SOFS8 2.8 umol/kg
%PO4 ave from RAS SOFS8 0.9
%NH4 set to 2 umol/kg
%H2S set to 0

fCO2_float=[];
for i=1:length(Alk_LIAR_20)
    [DATA, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    fCO2_float(end+1,i) = DATA(23); % uatm
end


% now combine this with atmospheric fCO2
%[F_CO2]=FCO2_CWE(DfCO2,T,S,u);