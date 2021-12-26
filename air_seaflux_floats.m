% convert float salinity into alkalinity and then airsea flux
% first go find the float data
% float IDs, choose one: 5906623 (SOTS), 5906624 (55S)
clear all
close all

float_ID = '5906623'

addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0\seawater_ver3_0'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper'

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)



% SOTS float adjustment parameters
adj_LIR_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_LIR_Deep');
adj_LIR_Shallow = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_LIR_Shallow');
adj_Williams_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_Williams_Deep');


% float_ID = '5906624' % 55S
% S55 float adjustment parameters
% adj = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55');



% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name]

% read necessary variables
lat = ncread(file,'LATITUDE');
lon = ncread(file,'LONGITUDE');
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
JULD = ncread(file,'JULD');
pres = ncread(file,'PRES'); %dbar
temp = ncread(file, 'TEMP_ADJUSTED');
temp_QC = ncread(file, 'TEMP_ADJUSTED_QC');
oxy = ncread(file,'DOXY_ADJUSTED');
oxy_QC = ncread(file,'DOXY_ADJUSTED_QC');
NO3 = ncread(file,'NITRATE');
%pH = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH = ncread(file,'PH_IN_SITU_TOTAL');
% pH_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
N = ncread(file, 'NITRATE_ADJUSTED');
S = ncread(file,'PSAL');
S_QC = ncread(file,'PSAL_QC');
Cycle = ncread(file,'CYCLE_NUMBER');
fs = size(pH);

SOTS_float_data = [];
SOTS_float_data.TEMP = temp;
SOTS_float_data.pH = pH;
SOTS_float_data.NO3 = NO3;
SOTS_float_data.pres = pres;
SOTS_float_data.psal = S;

% apply float correctsion as per Maurer, 2021 (the wrong, but quick way of
% adjusting pH rather than k0)

for n = 1:fs(2)
    idx_L_D = find(Cycle(n)>=adj_LIR_Deep.Cycle_start(:) & Cycle(n)<=adj_LIR_Deep.cycle_end(:));
    idx_L_S = find(Cycle(n)>=adj_LIR_Shallow.Cycle_start(:) & Cycle(n)<=adj_LIR_Shallow.cycle_end(:));
    idx_W_D = find(Cycle(n)>=adj_Williams_Deep.Cycle_start(:) & Cycle(n)<=adj_Williams_Deep.cycle_end(:));
    
    SOTS_float_data.pH_LIR_Deep(:,n) = pH(:,n)-adj_LIR_Deep.offset(idx_L_D)-adj_LIR_Deep.drift(idx_L_D)*(JULD(n)-JULD(adj_LIR_Deep.Cycle_start(idx_L_D)))/365;
    SOTS_float_data.pH_LIR_Shallow(:,n) = pH(:,n)-adj_LIR_Shallow.offset(idx_L_S)-adj_LIR_Shallow.drift(idx_L_S)*(JULD(n)-JULD(adj_LIR_Shallow.Cycle_start(idx_L_S)))/365;
    SOTS_float_data.pH_Williams_Deep(:,n) = pH(:,n)-adj_Williams_Deep.offset(idx_W_D)-adj_Williams_Deep.drift(idx_W_D)*(JULD(n)-JULD(adj_Williams_Deep.Cycle_start(idx_W_D)))/365;

end
% load the Cape Grim monthly xCO2 data (in dry air) and interpolate
% with a piecewise cubic Hermite spline (as per Gray paper)to hourly 
% data

file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\';

file = [file_path 'CapeGrim_CO2_data_download_18122021.xlsx'];

CGdata = readtable(file, 'Sheet', 'CapeGrim_CO2_data_download');
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
% file = [file_path '\S55_adaptor.mars.internal-1639715049.643984-2437-12-3c2e7bca-d3f8-4d73-a52d-ede241f56e62.nc'];
file = [file_path '\SOTS_float_adaptor.mars.internal-1635986988.5401623-6576-5-bb73848d-a029-46b4-abef-e1d94760fbeb.nc'];


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
Alk_ES = [];
Alk_ES_20 = [];
SOTS_float_data.Alk_pres=[];
pH_L_D_20 =[];
pH_L_S_20 =[];
pH_W_D_20 =[];
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

    Meas = [S(msk,i) oxy(msk,i) temp(msk,i)];
    
    % define the measurement IDs (salinity, oxygen, temperature)
    MeasID = [1 6 7];

    % run the LIAR code to create total alkalinity     
    Alk_LIAR(1:size(oxy(msk,i),1),end+1) = LIAR(Coor, Meas, MeasID);
    
    % calculate Alkalinity based on the linear regression for the top 100m
    % from the RAS alk and dic QC report
    % alk (umol/kg) = 39.23* salinity + 937.3
    Alk_ES(1:size(oxy(msk,i),1),end+1)  = 39.23*S(i)+ 937.3;
        
    % now we just want the top 20m average for the CO2SYS calcluations
    idx_L_D = pres(msk,i)<20;
    % a complicated way of getting an index vector for the accumarray
    % function
    c=ones(length(idx_L_D),1);
    c=c(idx_L_D);
    b=zeros(length(Alk_LIAR(:,1)),1);
    b(c~=0)=c;
    b=b+1;
    
    e=zeros(length(pres(msk,i)),1);
    e(c~=0)=c;
    e=e+1;
    
    % we only want the first value of accumarray (i.e. the first 20m)  
    a = accumarray(b,Alk_LIAR(:,i),[],@(x)mean(x,'omitnan'));
    Alk_LIAR_20(end+1,1) = a(1);
    a_ES = accumarray(b,Alk_ES(:,1),[],@(x)mean(x,'omitnan'));
    Alk_ES_20(end+1,1) = a_ES(1);
    p_L_D = accumarray(e,SOTS_float_data.pH_LIR_Deep(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_D_20(end+1,1) = p_L_D(1);
    p_L_S = accumarray(e,SOTS_float_data.pH_LIR_Shallow(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_S_20(end+1,1) = p_L_S(1);
    p_W_D = accumarray(e,SOTS_float_data.pH_Williams_Deep(msk,i),[],@(x)mean(x,'omitnan'));
    pH_W_D_20(end+1,1) = p_W_D(1);
    s = accumarray(e,S(msk,i),[],@(x)mean(x,'omitnan'));
    S_20(end+1,1) = s(1);
    pp = accumarray(e,pres(msk,i),[],@(x)mean(x,'omitnan'));
    pres_20(end+1,1) = pp(1);
    t = accumarray(e,temp(msk,i),[],@(x)mean(x,'omitnan'));
    temp_20(end+1,1) = t(1);
    SOTS_float_data.Alk_pres(1:size(oxy(msk,i),1),end+1) = pres(msk,i);
    
end
Alk_LIAR(Alk_LIAR==0)=NaN;
Alk_ES(Alk_ES==0)=NaN;

SOTS_float_data.Alk_LIAR = Alk_LIAR;
SOTS_float_data.Alk_ES = Alk_ES;


%Si ave from RAS SOFS8 2.8 umol/kg
%PO4 ave from RAS SOFS8 0.9
%NH4 set to 2 umol/kg
%H2S set to 0

% calculate fCO2 for float profiles with CO2SYS
%SOTS_float_data=[];
for i=1:length(Alk_LIAR_20)
    [DATA_L_D, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_L_D(i) = DATA_L_D(23); % uatm
    [DATA_L_S, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_S_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_L_S(i) = DATA_L_S(23); % uatm
    [DATA_W_D, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_W_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_W_D(i) = DATA_W_D(23); % uatm 
    
    [DATA_L_D_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_L_D_ES(i) = DATA_L_D_ES(23); % uatm
    [DATA_L_S_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_L_S_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_L_S_ES(i) = DATA_L_S_ES(23); % uatm
    [DATA_W_D_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_W_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    SOTS_float_data.fCO2_W_D_ES(i) = DATA_W_D_ES(23); % uatm 
    
    SOTS_float_data.lat(i) = lat(i);
    SOTS_float_data.lon(i) = lon(i);
    SOTS_float_data.time(i) = time_float(i);
    SOTS_float_data.TEMP_20(i) = temp_20(i);
    SOTS_float_data.PSAL_20(i) = S_20(i);
    SOTS_float_data.pH_L_D_20(i) = pH_L_D_20(i);
    SOTS_float_data.pH_L_S_20(i) = pH_L_S_20(i);
    SOTS_float_data.pH_W_D_20(i) = pH_W_D_20(i);
    SOTS_float_data.Alk_20(i) = Alk_LIAR_20(i);
    
end

%bias correction of float pH = -0.034529 * pH(25C) + 0.26709
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
    DfCO2_L_D = SOTS_float_data.fCO2_L_D(i)-pCO2_atm;
    DfCO2_L_S = SOTS_float_data.fCO2_L_S(i)-pCO2_atm;
    DfCO2_W_D = SOTS_float_data.fCO2_W_D(i)-pCO2_atm;
    [F_CO2_float_L_D]=FCO2_CWE(DfCO2_L_D,t_2_f,S_20(i),wsp_f);
    [F_CO2_float_L_S]=FCO2_CWE(DfCO2_L_S,t_2_f,S_20(i),wsp_f);
    [F_CO2_float_W_D]=FCO2_CWE(DfCO2_W_D,t_2_f,S_20(i),wsp_f);
    DfCO2_L_D_ES = SOTS_float_data.fCO2_L_D_ES(i)-pCO2_atm;
    DfCO2_L_S_ES = SOTS_float_data.fCO2_L_S_ES(i)-pCO2_atm;
    DfCO2_W_D_ES = SOTS_float_data.fCO2_W_D_ES(i)-pCO2_atm;
    [F_CO2_float_L_D_ES]=FCO2_CWE(DfCO2_L_D_ES,t_2_f,S_20(i),wsp_f);
    [F_CO2_float_L_S_ES]=FCO2_CWE(DfCO2_L_S_ES,t_2_f,S_20(i),wsp_f);
    [F_CO2_float_W_D_ES]=FCO2_CWE(DfCO2_W_D_ES,t_2_f,S_20(i),wsp_f);
    
    
    SOTS_float_data.flux_L_D(i) = F_CO2_float_L_D;
    SOTS_float_data.flux_L_S(i) = F_CO2_float_L_S;
    SOTS_float_data.flux_W_D(i) = F_CO2_float_W_D;
    SOTS_float_data.flux_L_D_ES(i) = F_CO2_float_L_D_ES;
    SOTS_float_data.flux_L_S_ES(i) = F_CO2_float_L_S_ES;
    SOTS_float_data.flux_W_D_ES(i) = F_CO2_float_W_D_ES;
    
    
    SOTS_float_data.wsp(i) = wsp_f;
    
end




clearvars -except SOTS_float_data

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('SOTS_float_data.mat')



% figure()
% geoscatter(fCO2_float.lat, fCO2_float.lon, datenum(fCO2_float.time)/4000, fCO2_float.flux,'^')
% hold on
% geoscatter(lat_SOLACE, lon_SOLACE, datenum(time_SOLACE)/6000, F_CO2_SOLACE,'.')
% colorbar

% figure()
% scatter(fCO2_float.time, fCO2_float.flux,'o')
% 
% figure()
% subplot(2,1,1)
% title('SOTS float and SOLACE UW data')
% yyaxis left
% plot(time_float,lat,'-b')
% hold on
% plot(time_SOLACE,lat_SOLACE,'--b')%,'MarkerSize',2)
% ylabel('Latitude')
% hold off
% yyaxis right
% plot(time_float,lon,'-r')
% hold on
% plot(time_SOLACE,lon_SOLACE,'--r')%,'MarkerSize',2)
% ylabel('Longitude')
% hold off
% xlabel('Time')
% 
% subplot(2,1,2)
% yyaxis left
% plot(fCO2_float.time, fCO2_float.flux,'ob','MarkerSize',4)
% hold on
% plot(time_SOLACE,F_CO2_SOLACE,'-r')
% ylabel('Air-sea CO2 flux mmol m^-2 d^-1')
% hold off
% xlabel('Time')
