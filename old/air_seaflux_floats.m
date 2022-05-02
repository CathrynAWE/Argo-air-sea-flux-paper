% convert float salinity into alkalinity and then airsea flux
% first go find the float data
% float IDs, choose one: 5906623 (SOTS), 5906624 (55S)
clear all
% close all

float_ID = '5906623' % SOTS float
% float_ID = '5906624' % 55S float
%%%%% change structure name!
SOTS_float_data = [];

addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0\seawater_ver3_0'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper'

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% float pH adjustments
if float_ID == '5906623' % SOTS
    % SOTS float adjustment parameters
    adj_LIR_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_LIR_Deep');
    adj_LIR_Shallow = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_LIR_Shallow');
    adj_Williams_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','SOTS_Williams_Deep');
elseif float_ID == '5906624' % 55S
    %S55 float adjustment parameters
    adj_LIR_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55_LIR_Deep');
    adj_LIR_Shallow = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55_LIR_Shallow');
    adj_Williams_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55_Williams_Deep');
end


% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name];

% read necessary variables from float profile
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

%%%% clear out the bad profiles and all their successive calculations
if float_ID == '5906623'
    pH(:,10:11) = NaN;
    pH(:,23) = NaN;
elseif float_ID == '5906624'
    pH(:,1:3) = NaN;
end



SOTS_float_data.TEMP = temp;
SOTS_float_data.pH = pH;
SOTS_float_data.oxy = oxy;
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

file = [file_path 'CapeGrim_CO2_data_download.xlsx'];

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


%%%%% before pCO2 is calculated from float pH and estimated Alk
%%%%% a bias correction needs to be made to bring float pH (akin to
%%%%% photospectormetrically measured pH) - based pCO2 calculation
%%%%% in line with measured pCO2 (or pCO2 estimated from DIC and Alk)
%bias correction of float = -0.034529 * pH(25C) + 0.26709 Williams, 2017


%%%%% set the reference depth (1500m as per Williams, 2017)
%%%%% I will use this for the LD and WD corrected profiles
ref_depth_D = 1500;
%%% set the tolerance to the reference depth
ref_depth_t_D = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i= 1:size(SOTS_float_data.pH,2)
    pres_ref_id_D = SOTS_float_data.pres(:,i)>=(ref_depth_D-ref_depth_t_D) & SOTS_float_data.pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
    pH_ref.pres_D(i) = mean(SOTS_float_data.pres(pres_ref_id_D,i),'omitnan');
    %%% find the pH at the reference depths
    pH_ref.pH_D(i) = mean(SOTS_float_data.pH_LIR_Deep(pres_ref_id_D,i),'omitnan');
    pH_ref.temp_D(i) = mean(SOTS_float_data.TEMP(pres_ref_id_D,i),'omitnan');
    pH_ref.psal_D(i) = mean(SOTS_float_data.psal(pres_ref_id_D,i),'omitnan');
    %%% convert the pH at the reference depths to what it would be at 25C
    [pH_25_D] = CO2SYS(2290,pH_ref.pH_D(i),1,3,pH_ref.psal_D(i),...
    pH_ref.temp_D(i),25, pH_ref.pres_D(i),0,2.8,0.9,2,0,1,10,1,2,2);
    pH_ref.pH_25_D(i) = pH_25_D(21);
    %%% take out the -999 and replace with NaN for subsequent calculations
    if pH_ref.pH_25_D(i) == -999
        pH_ref.pH_25_D(i) = NaN;
    end
%     %%% calculate the bias correction
    pH_ref.pH_corr_D(i) = -0.034529 * pH_ref.pH_25_D(i) + 0.26709;
end

%%% find the profile numbers where we have a deep enough measurement
pH_ref.profile_no_D = find(~isnan(pH_ref.pH_D));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% set the reference depth 950m this time for the LS profiles
%%%%% I will use this for the LS corrected profiles

ref_depth_S = 950;
%%% set the tolerance to the reference depth
ref_depth_t_S = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i= 1:size(SOTS_float_data.pH,2)
    pres_ref_id_S = SOTS_float_data.pres(:,i)>=(ref_depth_S-ref_depth_t_S) & SOTS_float_data.pres(:,i)<=(ref_depth_S+ref_depth_t_S); 
    pH_ref.pres_S(i) = mean(SOTS_float_data.pres(pres_ref_id_S,i),'omitnan');
    %%% find the pH at the reference depths
    pH_ref.pH_S(i) = mean(SOTS_float_data.pH_LIR_Deep(pres_ref_id_S,i),'omitnan');
    pH_ref.temp_S(i) = mean(SOTS_float_data.TEMP(pres_ref_id_S,i),'omitnan');
    pH_ref.psal_S(i) = mean(SOTS_float_data.psal(pres_ref_id_S,i),'omitnan');
    %%% convert the pH at the reference depths to what it would be at 25C
    [pH_25_S] = CO2SYS(2290,pH_ref.pH_S(i),1,3,pH_ref.psal_S(i),...
    pH_ref.temp_S(i),25, pH_ref.pres_S(i),0,2.8,0.9,2,0,1,10,1,2,2);
    pH_ref.pH_25_S(i) = pH_25_S(21);
    %%% take out the -999 and replace with NaN for subsequent calculations
    if pH_ref.pH_25_S(i) == -999
        pH_ref.pH_25_S(i) = NaN;
    end
%     %%% calculate the bias correction
    pH_ref.pH_corr_S(i) = -0.034529 * pH_ref.pH_25_S(i) + 0.26709;
end

%%% find the profile numbers where we have a deep enough measurement
pH_ref.profile_no_S = find(~isnan(pH_ref.pH_S));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% now each pH profile has to be corrected by adding the bias
%%% I will use the correction that is closest to each profile
for i = 1:size(SOTS_float_data.pH,2)
%for i=15   
    ind = knnsearch(pH_ref.profile_no_D(:),i);
    idx = pH_ref.profile_no_D(ind);
    SOTS_float_data.pH_LD_corr(:,i) = SOTS_float_data.pH_LIR_Deep(:,i)+pH_ref.pH_corr_D(idx);
    SOTS_float_data.pH_WD_corr(:,i) = SOTS_float_data.pH_Williams_Deep(:,i)+pH_ref.pH_corr_D(idx);
end

for i = 1:size(SOTS_float_data.pH,2)
    ind = knnsearch(pH_ref.profile_no_S(:),i);
    idx = pH_ref.profile_no_S(ind);
    SOTS_float_data.pH_LS_corr(:,i) = SOTS_float_data.pH_LIR_Shallow(:,i)+pH_ref.pH_corr_S(idx);
end

% convert float Salinity, oxygen and temperature into alkalinity with LIAR
Alk_LIAR=[];
Alk_LIAR_20=[];
Alk_ES = [];
Alk_ES_20 = [];
SOTS_float_data.Alk_pres=[];
pH_L_D_20 =[];
pH_L_S_20 =[];
pH_W_D_20 =[];
pH_L_D_corr_20 =[];
pH_L_S_corr_20 =[];
pH_W_D_corr_20 =[];

S_20 = [];
temp_20 = [];
pres_20 =[];

for i = 1:fs(2)
%for i =1:3
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
    Alk_ES(1:size(oxy(msk,i),1),end+1)  = 39.23*S(msk,i)+ 937.3;
        
    % now we just want the top 20m average for the CO2SYS calcluations
    idx = pres(msk,i)<=20;
    % get an index vector for the accumarray function
    c=ones(length(idx),1)+idx;
    len=length(Alk_LIAR(:,i))-length(c);
    b=wextend('1','zpd',c,len,'r');
    b(b==0)=1;
    
    % we only want the first value of accumarray (i.e. the first 20m)  
    a = accumarray(b,Alk_LIAR(:,i),[],@(x)mean(x,'omitnan'));
    Alk_LIAR_20(end+1,1) = a(2);
    a_ES = accumarray(b,Alk_ES(:,1),[],@(x)mean(x,'omitnan'));
    Alk_ES_20(end+1,1) = a_ES(2);
    p_L_D = accumarray(c,SOTS_float_data.pH_LIR_Deep(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_D_20(end+1,1) = p_L_D(2);
    p_L_S = accumarray(c,SOTS_float_data.pH_LIR_Shallow(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_S_20(end+1,1) = p_L_S(2);
    p_W_D = accumarray(c,SOTS_float_data.pH_Williams_Deep(msk,i),[],@(x)mean(x,'omitnan'));
    pH_W_D_20(end+1,1) = p_W_D(2);
    p_L_D_corr = accumarray(c,SOTS_float_data.pH_LD_corr(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_D_corr_20(end+1,1) = p_L_D_corr(2); 
    p_L_S_corr = accumarray(c,SOTS_float_data.pH_LS_corr(msk,i),[],@(x)mean(x,'omitnan'));
    pH_L_S_corr_20(end+1,1) = p_L_S_corr(2); 
    p_W_D_corr = accumarray(c,SOTS_float_data.pH_WD_corr(msk,i),[],@(x)mean(x,'omitnan'));
    pH_W_D_corr_20(end+1,1) = p_W_D_corr(2); 
    
    s = accumarray(c,S(msk,i),[],@(x)mean(x,'omitnan'));
    S_20(end+1,1) = s(2);
    pp = accumarray(c,pres(msk,i),[],@(x)mean(x,'omitnan'));
    pres_20(end+1,1) = pp(2);
    t = accumarray(c,temp(msk,i),[],@(x)mean(x,'omitnan'));
    temp_20(end+1,1) = t(2);
    SOTS_float_data.Alk_pres(1:size(oxy(msk,i),1),end+1) = pres(msk,i);
    
end
% there are trailing 0s in the alkalinity results because the vectors are
% not equally long
Alk_LIAR(Alk_LIAR==0)=NaN;
Alk_ES(Alk_ES==0)=NaN;

SOTS_float_data.Alk_LIAR = Alk_LIAR;
SOTS_float_data.Alk_ES = Alk_ES;



%Si ave from RAS SOFS8 2.8 umol/kg - WOA nuts, 2.5 for SOTS, 6.9 for S55
%PO4 ave from RAS SOFS8 0.9 - WOA nuts 0.9 for SOTS, 1.4 for S55
%NH4 set to 2 umol/kg
%H2S set to 0

% calculate fCO2 for float profiles with CO2SYS
%SOTS_float_data=[];
% I will use pCO2 for the floats, since I have pCO2 from Cape Grim
% I will use annual mean nutrients from WOA for the lat and lon range that
% the float covered in the year 2021

if float_ID == '5906623' % SOTS float
    load('WOA_nut_SOTS.mat')
    nut_data = WOA_nut_SOTS;
elseif float_ID == '5906624' % 55S float
    load('WOA_nut_S55.mat')
    nut_data = WOA_nut_S55;
end
    
for i=1:length(Alk_LIAR_20)
    [DATA_L_D, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_D(i) = DATA_L_D(22); % uatm
    
    [DATA_L_D_20, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),20,20,nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_D_20m(i) = DATA_L_D_20(22); % uatm
    
    [DATA_L_S, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_S_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_S(i) = DATA_L_S(22); % uatm
    [DATA_W_D, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_W_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_W_D(i) = DATA_W_D(22); % uatm 
    
    [DATA_L_D_corr, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_corr_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_D_corr(i) = DATA_L_D_corr(22); % uatm
    [DATA_L_S_corr, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_S_corr_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_S_corr(i) = DATA_L_S_corr(22); % uatm
    [DATA_W_D_corr, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_W_D_corr_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_W_D_corr(i) = DATA_W_D_corr(22); % uatm

    [DATA_L_D_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_D_ES(i) = DATA_L_D_ES(22); % uatm
    [DATA_L_S_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_L_S_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_L_S_ES(i) = DATA_L_S_ES(22); % uatm
    [DATA_W_D_ES, HEADERS, NICEHEADERS]  = CO2SYS(Alk_ES_20(i),pH_W_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    SOTS_float_data.pCO2_W_D_ES(i) = DATA_W_D_ES(22); % uatm 
    
    SOTS_float_data.lat(i) = lat(i);
    SOTS_float_data.lon(i) = lon(i);
    SOTS_float_data.time(i) = time_float(i);
    SOTS_float_data.TEMP_20(i) = temp_20(i);
    SOTS_float_data.PSAL_20(i) = S_20(i);
    SOTS_float_data.pH_L_D_20(i) = pH_L_D_20(i);
    SOTS_float_data.pH_L_S_20(i) = pH_L_S_20(i);
    SOTS_float_data.pH_W_D_20(i) = pH_W_D_20(i);
    SOTS_float_data.pH_L_D_20_corr(i) = pH_L_D_corr_20(i);
    SOTS_float_data.pH_L_S_20_corr(i) = pH_L_S_corr_20(i);
    SOTS_float_data.pH_W_D_20_corr(i) = pH_W_D_corr_20(i);
        
    SOTS_float_data.Alk_LIAR_20(i) = Alk_LIAR_20(i);
    SOTS_float_data.Alk_ES_20(i) = Alk_ES_20(i);
    SOTS_float_data.pres_20(i) = pres_20(i);
end

SOTS_float_data.pCO2_L_D(SOTS_float_data.pCO2_L_D== -999) = NaN;
SOTS_float_data.pCO2_L_D_20m(SOTS_float_data.pCO2_L_D_20m== -999) = NaN;
SOTS_float_data.pCO2_L_S(SOTS_float_data.pCO2_L_S == -999) = NaN;
SOTS_float_data.pCO2_W_D(SOTS_float_data.pCO2_W_D== -999) = NaN;
SOTS_float_data.pCO2_L_D_corr(SOTS_float_data.pCO2_L_D_corr== -999) = NaN;
SOTS_float_data.pCO2_L_S_corr(SOTS_float_data.pCO2_L_S_corr== -999) = NaN;
SOTS_float_data.pCO2_W_D_corr(SOTS_float_data.pCO2_W_D_corr == -999) = NaN;
SOTS_float_data.pCO2_L_D_ES(SOTS_float_data.pCO2_L_D_ES== -999) = NaN;
SOTS_float_data.pCO2_L_S_ES(SOTS_float_data.pCO2_L_S_ES== -999) = NaN;
SOTS_float_data.pCO2_W_D_ES(SOTS_float_data.pCO2_W_D_ES== -999) = NaN;


% read the ERA5 data for the area of the float
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript';
if float_ID =='5906623'
    file = [file_path '\SOTS_float_adaptor.mars.internal-1635986988.5401623-6576-5-bb73848d-a029-46b4-abef-e1d94760fbeb.nc'];
elseif float_ID =='5906624'
    file = [file_path '\S55_adaptor.mars.internal-1639715049.643984-2437-12-3c2e7bca-d3f8-4d73-a52d-ede241f56e62.nc'];
end

u_10 = ncread(file,'u10');
v_10 = ncread(file, 'v10');
t_2 = ncread(file, 't2m')-273.15; % Kelvin converted to Celsius
msp = ncread(file, 'msl');
EraTime = hours(ncread(file, 'time'))+ datetime(1900,1,1);
EraLat = ncread(file, 'latitude');
EraLon = ncread(file, 'longitude');

if float_ID == '5906623'
    load('NCEP_SOTS_float.mat')
elseif float_ID == '5906624'
    load ('NCEP_S55_float.mat')
end

disp('calculating windspeeds')

for i = 1:fs(2)
%for i = 156   
    % find the ERA file index that is closest in time and space to the float
    idx_T = knnsearch(datenum(EraTime(:)),datenum(time_float(i)));
    delta_time = (EraTime(idx_T)-time_float(i))/24;
    idx_lat = knnsearch(EraLat(:),lat(i));
    idx_lon = knnsearch(EraLon(:),lon(i));
    
        
    wsp_f = sqrt(u_10(idx_lon,idx_lat,1,idx_T)^2 + v_10(idx_lon,idx_lat,1,idx_T)^2);
    %m s-2
    t_2_f = t_2(idx_lon,idx_lat,1,idx_T); %C
    msp_f = msp(idx_lon,idx_lat,1,idx_T); %Pa
    if abs(delta_time)>1
        wsp_f = NaN;
        t_2_f = NaN;
        msp_f = NaN;
    end
        
    if isnan(wsp_f)
        wsp_f = float_ncep.wnd_ncep_f(i);
    end
    
    if isnan(t_2_f)
        t_2_f = float_ncep.airT_ncep_f(i);
    end
        
    if isnan(msp_f)
        msp_f = float_ncep.mslp_ncep_f(i);
    end
        
    
    % water vapour pressure 6.11 * 10^((7.5*(T))/(237.3+(T)))
    % Temp in degrees C 
    pH20_f = 6.11 * 10^((7.5*t_2_f)/(273.3+t_2_f)); %mbar = hPa

    % convert xCO2 in dry air to pCO2
    idx_C = knnsearch(datenum(CGdata_interp.spl_time(:)),datenum(time_float(i)));
    %uatm
    pCO2_uatm = CGdata_interp.ppm_spl(idx_C) * ((msp_f/100 - pH20_f))*(9.8692326671601*10^-4); %uatm
    
    % now combine this with atmospheric fCO2
    DpCO2_L_D = SOTS_float_data.pCO2_L_D(i)-pCO2_uatm;
    DpCO2_L_D_20m = SOTS_float_data.pCO2_L_D_20m(i)-pCO2_uatm;
    DpCO2_L_S = SOTS_float_data.pCO2_L_S(i)-pCO2_uatm;
    DpCO2_W_D = SOTS_float_data.pCO2_W_D(i)-pCO2_uatm;
    DpCO2_L_D_corr = SOTS_float_data.pCO2_L_D_corr(i)-pCO2_uatm;
    DpCO2_L_S_corr = SOTS_float_data.pCO2_L_S_corr(i)-pCO2_uatm;
    DpCO2_W_D_corr = SOTS_float_data.pCO2_W_D_corr(i)-pCO2_uatm;
    
    [F_CO2_float_L_D]=FCO2_CWE(DpCO2_L_D,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_D_20m]=FCO2_CWE(DpCO2_L_D_20m,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_S]=FCO2_CWE(DpCO2_L_S,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_W_D]=FCO2_CWE(DpCO2_W_D,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_D_corr]=FCO2_CWE(DpCO2_L_D_corr,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_S_corr]=FCO2_CWE(DpCO2_L_S_corr,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_W_D_corr]=FCO2_CWE(DpCO2_W_D_corr,temp_20(i),S_20(i),wsp_f);


%     [F_CO2_float_W_D]=FCO2_CWE(DpCO2_W_D,t_2_f,S_20(i),wsp_f);
    DpCO2_L_D_ES = SOTS_float_data.pCO2_L_D_ES(i)-pCO2_uatm;
    DpCO2_L_S_ES = SOTS_float_data.pCO2_L_S_ES(i)-pCO2_uatm;
    DpCO2_W_D_ES = SOTS_float_data.pCO2_W_D_ES(i)-pCO2_uatm;
%     [F_CO2_float_L_D_ES_airtemp]=FCO2_CWE(DpCO2_L_D_ES,t_2_f,S_20(i),wsp_f);
    [F_CO2_float_L_D_ES]=FCO2_CWE(DpCO2_L_D_ES,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_S_ES]=FCO2_CWE(DpCO2_L_S_ES,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_W_D_ES]=FCO2_CWE(DpCO2_W_D_ES,temp_20(i),S_20(i),wsp_f);
    
    
    SOTS_float_data.flux_L_D(i) = F_CO2_float_L_D;
    SOTS_float_data.flux_L_D_20m(i) = F_CO2_float_L_D_20m;
    SOTS_float_data.flux_L_S(i) = F_CO2_float_L_S;
    SOTS_float_data.flux_W_D(i) = F_CO2_float_W_D;
    SOTS_float_data.flux_L_D_corr(i) = F_CO2_float_L_D_corr;
    SOTS_float_data.flux_L_S_corr(i) = F_CO2_float_L_S_corr;
    SOTS_float_data.flux_W_D_corr(i) = F_CO2_float_W_D_corr;
    

    SOTS_float_data.flux_L_D_ES(i) = F_CO2_float_L_D_ES;
%     SOTS_float_data.flux_L_D_ES_airtemp(i) = F_CO2_float_L_D_ES_airtemp;
    SOTS_float_data.flux_L_S_ES(i) = F_CO2_float_L_S_ES;
    SOTS_float_data.flux_W_D_ES(i) = F_CO2_float_W_D_ES;
    
    SOTS_float_data.pCO2_uatm(i) = pCO2_uatm;
    SOTS_float_data.CG_xCO2_uatm(i) = CGdata_interp.ppm_spl(idx_C);
    SOTS_float_data.wsp(i) = wsp_f;
    SOTS_float_data.t_2(i) = t_2_f;
    SOTS_float_data.msp(i) = msp_f;
end

%%%% now I just want 2021, because it is a full year
SOTS_float_data.Month = month(SOTS_float_data.time);
SOTS_float_data.year = year(SOTS_float_data.time);

SOTS_float_data.data_comp_2021 = table(SOTS_float_data.Month(SOTS_float_data.year==2021)',...
    SOTS_float_data.time(SOTS_float_data.year==2021)',...
    SOTS_float_data.pH_L_D_20(SOTS_float_data.year==2021)',SOTS_float_data.pH_L_S_20(SOTS_float_data.year==2021)',...
    SOTS_float_data.pH_W_D_20(SOTS_float_data.year==2021)',SOTS_float_data.pH_L_D_20_corr(SOTS_float_data.year==2021)',...
    SOTS_float_data.pH_L_S_20_corr(SOTS_float_data.year==2021)',SOTS_float_data.pH_W_D_20_corr(SOTS_float_data.year==2021)',...
    SOTS_float_data.Alk_LIAR_20(SOTS_float_data.year==2021)',SOTS_float_data.Alk_ES_20(SOTS_float_data.year==2021)',...
    SOTS_float_data.flux_L_D(SOTS_float_data.year==2021)', SOTS_float_data.flux_L_S(SOTS_float_data.year==2021)',...
    SOTS_float_data.flux_W_D(SOTS_float_data.year==2021)',SOTS_float_data.flux_L_D_corr(SOTS_float_data.year==2021)',...
    SOTS_float_data.flux_L_S_corr(SOTS_float_data.year==2021)',SOTS_float_data.flux_W_D_corr(SOTS_float_data.year==2021)',...
    SOTS_float_data.flux_L_D_ES(SOTS_float_data.year==2021)', SOTS_float_data.flux_L_S_ES(SOTS_float_data.year==2021)',...
    SOTS_float_data.flux_W_D_ES(SOTS_float_data.year==2021)', 'VariableNames',{'Month','Time','pH_LD_20','pH_LS_20',...
    'pH_WD_20','pH_LD_20_corr','pH_LS_20_corr','pH_WD_20_corr','Alk_LIAR_20','Alk_ES_20',...
    'flux_LD','flux_LS','flux_WD','flux_LD_corr','flux_LS_corr','flux_WD_corr'...
    'flux_LD_ES','flux_LS_ES','flux_WD_ES'});

SOTS_float_data.pH_LD_20_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_LD_20,[],@(x)mean(x,'omitnan'));
SOTS_float_data.pH_LS_20_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_LS_20,[],@(x)mean(x,'omitnan'));
SOTS_float_data.pH_WD_20_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_WD_20,[],@(x)mean(x,'omitnan'));
SOTS_float_data.pH_LD_20_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_LD_20_corr,[],@(x)mean(x,'omitnan'));
SOTS_float_data.pH_LS_20_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_LS_20_corr,[],@(x)mean(x,'omitnan'));
SOTS_float_data.pH_WD_20_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.pH_WD_20_corr,[],@(x)mean(x,'omitnan'));

SOTS_float_data.Alk_LIAR_20_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.Alk_LIAR_20,[],@(x)mean(x,'omitnan'));
SOTS_float_data.Alk_ES_20_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.Alk_ES_20,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_LD_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LD,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_LS_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LS,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_WD_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_WD,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_LD_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LD_corr,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_LS_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LS_corr,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_WD_corr_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_WD_corr,[],@(x)mean(x,'omitnan'));

SOTS_float_data.flux_LD_ES_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LD_ES,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_LS_ES_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LS_ES,[],@(x)mean(x,'omitnan'));
SOTS_float_data.flux_WD_ES_mo_ave = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_WD_ES,[],@(x)mean(x,'omitnan'));

SOTS_float_data.flux_LD_corr_mo_SD = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LD_corr,[],@(x)std(x,'omitnan'));
SOTS_float_data.flux_LS_corr_mo_SD = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_LS_corr,[],@(x)std(x,'omitnan'));
SOTS_float_data.flux_WD_corr_mo_SD = accumarray(SOTS_float_data.data_comp_2021.Month,SOTS_float_data.data_comp_2021.flux_WD_corr,[],@(x)std(x,'omitnan'));


SOTS_float_data.mo_ave_month = unique(SOTS_float_data.data_comp_2021.Month);

% for the SOTS float only
% calculate the distance between float and mooring

if float_ID == '5906623'
    load('mooring_data.mat')

    for i = 1:fs(2)
    %for i = 156   
        % find the mooring_data index that is closest in time to the float time
        idx_T = knnsearch(datenum(mooring_data.xCO2_time(:)),datenum(time_float(i)));
        idx_lat = knnsearch(mooring_data.xCO2_lat(:),lat(i));
        idx_lon = knnsearch(mooring_data.xCO2_lon(:),lon(i));

        SOTS_float_data.mooring_lat(i) = mooring_data.xCO2_lat(idx_lat);
        SOTS_float_data.mooring_lon(i) = mooring_data.xCO2_lon(idx_lon);

    end
end

if float_ID == '5906623'
    clearvars -except SOTS_float_data pH_ref CGdata_interp float_ID
elseif float_ID == '5906624'
    clearvars -except S55_float_data pH_ref CGdata_interp float_ID
end


path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
if float_ID == '5906623'
    save('SOTS_float_data.mat')
elseif float_ID == '5906624'
    save('S55_float_data.mat')
end


% figure()
% geoscatter(fCO2_float.lat, fCO2_float.lon, datenum(fCO2_float.time)/4000, fCO2_float.flux,'^')
% hold on
% geoscatter(lat_SOLACE, lon_SOLACE, datenum(time_SOLACE)/6000, F_CO2_SOLACE,'.')
% colorbar