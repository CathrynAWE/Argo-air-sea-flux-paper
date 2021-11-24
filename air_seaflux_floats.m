% convert float salinity into alkalinity and then airsea flux
% first go find the float data
% float IDs, choose one: 5906623 (SOTS), 5906624 (55S)
%clear all

float_ID = '5906623'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0\seawater_ver3_0'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper'

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name]

% read necessary variables
lat = ncread(file,'LATITUDE');
lon = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
pres = ncread(file,'PRES');
temp = ncread(file, 'TEMP_ADJUSTED');
oxy = ncread(file,'DOXY_ADJUSTED');
NO3 = ncread(file,'NITRATE');
%pH = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH = ncread(file,'PH_IN_SITU_TOTAL');
pH_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
N = ncread(file, 'NITRATE_ADJUSTED');
T = ncread(file,'TEMP');
S = ncread(file,'PSAL');
fs = size(pH);

% load the Cape Grim monthly xCO2 data (in dry air) and interpolate to
% hourly data
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\';

file = [file_path 'CapeGrim_CO2_data_download_08112021.xlsx'];

CGdata = readtable(file, 'Sheet', 'CapeGrim_CO2_data_download_0811');
% clean out the NAN lines (which were text in the original)
CGdata(isnan(CGdata.YYYY),:)=[];

CGdata.datestr = strcat(num2str(CGdata.YYYY),{'-'},num2str(CGdata.MM),{'-'}, num2str(CGdata.DD));

date = datetime(CGdata.datestr, 'InputFormat','yyyy-MM-dd');

CGdata_interp = [];
CGdata_interp.date = date;
CGdata_interp.ppm = CGdata.CO2_ppm_;
CGdata_interp.datehours=(datenum(date)-datenum('1990-1-1'))*24;
xq2 = CGdata_interp.datehours(1):1:CGdata_interp.datehours(end);
CGdata_interp.ppm_spl = pchip(CGdata_interp.datehours,CGdata.CO2_ppm_,xq2);
CGdata_interp.spl_time = datetime((xq2/24)+datenum('1990-1-1'),'ConvertFrom','datenum');

% read the ERA5 data for the area of the float
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript'
file = [file_path '\adaptor.mars.internal-1635986988.5401623-6576-5-bb73848d-a029-46b4-abef-e1d94760fbeb.nc'];

u_10 = ncread(file,'u10');
v_10 = ncread(file, 'v10');
t_2 = ncread(file, 't2m')-273.16; % Kelvin converted to Celsius
msp = ncread(file, 'msl');
EraTime = hours(ncread(file, 'time'))+ datetime(1900,1,1);
EraLat = ncread(file, 'latitude');
EraLon = ncread(file, 'longitude');



% convert float Salinity, oxygen and temperature into alkalinity with LIAR
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

% calculate fCO2 for float profiles with CO2SYS
fCO2_float=[];
for i=1:length(Alk_LIAR_20)
    [DATA, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    fCO2_float.fCO2(i) = DATA(23); % uatm
    fCO2_float.lat(i) = lat(i);
    fCO2_float.lon(i) = lon(i);
    fCO2_float.time(i) = time_float(i);
end



for i = 1:fs(2)
    
    % find the ERA file index that is closest in time to the float time
    idx_T = knnsearch(datenum(EraTime(:)),datenum(time_float(i)));
    idx_lat = knnsearch(EraLat(:),lat(i));
    idx_lon = knnsearch(EraLon(:),lon(i));

    wsp_f = sqrt(u_10(idx_lon,idx_lat,1,idx_T)^2 + v_10(idx_lon,idx_lat,1,idx_T)^2);
    t_2_f = t_2(idx_lon,idx_lat,1,idx_T);
    msp_f = msp(idx_lon,idx_lat,1,idx_T);

    % water vapour pressure 6.11 * 10^((7.5*(T-273.16))/(237.3+(T-273.16)))
    pH20_f = 6.11 * 10^((7.5*t_2_f)/(273.3+t_2_f));

    % convert xCO2 in dry air to pCO2
    idx_C = knnsearch(datenum(CGdata_interp.spl_time(:)),datenum(time_float(i)));
    %uatm
    pCO2_atm = CGdata_interp.ppm_spl(idx_C) * ((msp_f - pH20_f*100))*(9.8692326671601*10^-6);

    % now combine this with atmospheric fCO2
    DfCO2 = fCO2_float.fCO2(i)-pCO2_atm;

    [F_CO2]=FCO2_CWE(DfCO2,t_2_f,S_20(i),wsp_f);

    fCO2_float.flux(i) = F_CO2;
    
end


% figure()
% geoscatter(fCO2_float.lat, fCO2_float.lon, datenum(fCO2_float.time)/4000, fCO2_float.flux,'^')
% hold on
% geoscatter(lat_SOLACE, lon_SOLACE, datenum(time_SOLACE)/6000, F_CO2_SOLACE,'.')
% colorbar

% figure()
% scatter(fCO2_float.time, fCO2_float.flux,'o')

figure()
subplot(2,1,1)
title('SOTS float and SOLACE UW data')
yyaxis left
plot(time_float,lat,'-b')
hold on
plot(time_SOLACE,lat_SOLACE,'--b')%,'MarkerSize',2)
ylabel('Latitude')
hold off
yyaxis right
plot(time_float,lon,'-r')
hold on
plot(time_SOLACE,lon_SOLACE,'--r')%,'MarkerSize',2)
ylabel('Longitude')
hold off
xlabel('Time')

subplot(2,1,2)
yyaxis left
plot(fCO2_float.time, fCO2_float.flux,'ob','MarkerSize',4)
hold on
plot(time_SOLACE,F_CO2_SOLACE,'-r')
ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
hold off
xlabel('Time')
