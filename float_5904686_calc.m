%%%% this script is calculating corrected pH data for a float that was in
%%%% the vicinity to our current 55S float. We are trying to see how much
%%%% difference the various choices for correcting pH makes to say air sea
%%%% flux calculations. We also want to see how our corrections compare to
%%%% an expert's correction (Maurer's pH corrected in this file).
%%%% this float only has good data to cycle 101!

clear all
close all

float_ID = '5904686' % 55S complimentary float

addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0\seawater_ver3_0'
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper'

search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% float adjustment parameters
adj_LIR_Deep = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55_comp_LIR_Deep');
adj_LIR_Shallow = readtable('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\Float_adjustments.xlsx','Sheet','S55_comp_LIR_Shallow');

% find the Sprof
nc_files = dir('*.nc');
file = [nc_files.folder '\' nc_files.name];

% read necessary variables from float profile
lat = ncread(file,'LATITUDE');
lat = lat(1:101);
lon = ncread(file,'LONGITUDE');
lon =  lon(1:101);
time_float = ncread(file,'JULD')+ datetime(1950,1,1);
time_float = time_float(1:101);
JULD = ncread(file,'JULD');
JULD = JULD(1:101);
pres = ncread(file,'PRES'); %dbar
pres = pres(:,1:101);
temp = ncread(file, 'TEMP_ADJUSTED');
temp = temp(:,1:101);
temp_QC = ncread(file, 'TEMP_ADJUSTED_QC');
temp_QC = temp_QC(:,1:101);
oxy = ncread(file,'DOXY_ADJUSTED');
oxy = oxy(:,1:101);
oxy_QC = ncread(file,'DOXY_ADJUSTED_QC');
oxy_QC = oxy_QC(:,1:101);
pH_adj = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH_adj = pH_adj(:,1:101);
pH_adj_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
pH_adj_error = pH_adj_error(:,1:101);
pH_adj_QC = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_QC');
pH_adj_QC = pH_adj_QC(:,1:101);
pH_raw = ncread(file,'PH_IN_SITU_TOTAL');
pH_raw = pH_raw(:,1:101);
NO3 = ncread(file, 'NITRATE_ADJUSTED');
NO3 = NO3(:,1:101);
S = ncread(file,'PSAL');
S = S(:,1:101);
S_QC = ncread(file,'PSAL_QC');
S_QC = S_QC(:,1:101);
Cycle = ncread(file,'CYCLE_NUMBER');
Cycle = Cycle(1:101);
fs = size(pH_raw);

S55_comp_float_data = [];
S55_comp_float_data.TEMP = temp;
S55_comp_float_data.pH_raw = pH_raw;
S55_comp_float_data.pH_adj = pH_adj;
S55_comp_float_data.pH_adj_QC = pH_adj_QC;
S55_comp_float_data.pH_adj_error = pH_adj_error;
S55_comp_float_data.pres = pres;
S55_comp_float_data.psal = S;
S55_comp_float_data.psal_QC = S_QC;

% apply float correction as per Maurer, 2021 (the wrong, but quick way of
% adjusting pH rather than k0)
for n = 1:fs(2)
    idx_L_D = find(Cycle(n)>=adj_LIR_Deep.Cycle_start(:) & Cycle(n)<=adj_LIR_Deep.cycle_end(:));
    idx_L_S = find(Cycle(n)>=adj_LIR_Shallow.Cycle_start(:) & Cycle(n)<=adj_LIR_Shallow.cycle_end(:));
    
    S55_comp_float_data.pH_LIR_Deep(:,n) = pH_raw(:,n)-adj_LIR_Deep.offset(idx_L_D)-adj_LIR_Deep.drift(idx_L_D)*(JULD(n)-JULD(adj_LIR_Deep.Cycle_start(idx_L_D)))/365;
    S55_comp_float_data.pH_LIR_Shallow(:,n) = pH_raw(:,n)-adj_LIR_Shallow.offset(idx_L_S)-adj_LIR_Shallow.drift(idx_L_S)*(JULD(n)-JULD(adj_LIR_Shallow.Cycle_start(idx_L_S)))/365;
    
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
ref_depth = ones(1,size(S55_comp_float_data.pH_raw,2))*1500;
%%% set the tolerance to the reference depth
ref_depth_t = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i=1:size(S55_comp_float_data.pH_raw,2)
    pres_ref_id(i) = knnsearch(S55_comp_float_data.pres(:,i), ref_depth(i));
    pres_ref(i) = S55_comp_float_data.pres(pres_ref_id(i),i);
end

%%% take out the depth that are not within the depth tolerance
pres_ref(pres_ref >= (ref_depth(1) + ref_depth_t))=NaN;
pres_ref(pres_ref <= (ref_depth(1) - ref_depth_t))=NaN;

%%% find the profile numbers where we have a deep enough measurement
pH_ref_profile_ind=find(~isnan(pres_ref));

pH_ref_depth_ind=pres_ref_id(pH_ref_profile_ind);

%%% find the pH at the reference depths
for i=1:length(pH_ref_profile_ind)
    pH_ref.pH(i) = S55_comp_float_data.pH_raw(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.temp(i) = S55_comp_float_data.TEMP(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.pres(i) = S55_comp_float_data.pres(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.psal(i) = S55_comp_float_data.psal(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.profile_no(i) = pH_ref_profile_ind(i);
    %%%% calculate the pH at 25C and 0dbar
    [pH_25] = CO2SYS(2290,pH_ref.pH(i),1,3,pH_ref.psal(i),...
        pH_ref.temp(i),25, pH_ref.pres(i),0,2.8,0.9,2,0,1,10,1,2,2);
    pH_ref.pH_25(i) = pH_25(21);
    %%% calculate the bias correction
    pH_ref.pH_corr(i) = -0.034529 * pH_ref.pH_25(i) + 0.26709;
end

%%%% now each pH profile has to be corrected by adding the bias
%%% I will use the correction that is closest to each profile
for i = 1:size(S55_comp_float_data.pH_raw,2)
%for i=15   
    ind = knnsearch(pH_ref.profile_no(:),i);
    S55_comp_float_data.pH_LD_corr(:,i) = S55_comp_float_data.pH_LIR_Deep(:,i)+pH_ref.pH_corr(ind);
    S55_comp_float_data.pH_LS_corr(:,i) = S55_comp_float_data.pH_LIR_Shallow(:,i)+pH_ref.pH_corr(ind);
   
end


% convert float Salinity, oxygen and temperature into alkalinity with LIAR
Alk_LIAR=[];
Alk_LIAR_20=[];
S55_comp_float_data.Alk_pres=[];
pH_L_D_20 =[];
pH_L_S_20 =[];
pH_L_D_corr_20 =[];
pH_L_S_corr_20 =[];
pH_adj_20 = [];
S_20 = [];
temp_20 = [];
pres_20 =[];

for i = 1:fs(2)
%for i =1:3
    % mask out the NaNs in the float data, for the Alkalinity script
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
    if ~isempty(Meas)
        Alk_LIAR(1:size(oxy(msk,i),1),end+1) = LIAR(Coor, Meas, MeasID);
      
         
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
        p_L_D = accumarray(c,S55_comp_float_data.pH_LIR_Deep(msk,i),[],@(x)mean(x,'omitnan'));
        pH_L_D_20(end+1,1) = p_L_D(2);
        p_L_S = accumarray(c,S55_comp_float_data.pH_LIR_Shallow(msk,i),[],@(x)mean(x,'omitnan'));
        pH_L_S_20(end+1,1) = p_L_S(2);
        p_L_D_corr = accumarray(c,S55_comp_float_data.pH_LD_corr(msk,i),[],@(x)mean(x,'omitnan'));
        pH_L_D_corr_20(end+1,1) = p_L_D_corr(2); 
        p_L_S_corr = accumarray(c,S55_comp_float_data.pH_LS_corr(msk,i),[],@(x)mean(x,'omitnan'));
        pH_L_S_corr_20(end+1,1) = p_L_S_corr(2);
        p_adj_20 = accumarray(c,S55_comp_float_data.pH_adj(msk,i),[],@(x)mean(x,'omitnan'));
        pH_adj_20(end+1,1) = p_adj_20(2);
        

        s = accumarray(c,S(msk,i),[],@(x)mean(x,'omitnan'));
        S_20(end+1,1) = s(2);
        pp = accumarray(c,pres(msk,i),[],@(x)mean(x,'omitnan'));
        pres_20(end+1,1) = pp(2);
        t = accumarray(c,temp(msk,i),[],@(x)mean(x,'omitnan'));
        temp_20(end+1,1) = t(2);
        S55_comp_float_data.Alk_pres(1:size(oxy(msk,i),1),end+1) = pres(msk,i);
    else
        Alk_LIAR(1:153,end+1) = NaN;
        Alk_LIAR_20(end+1,1) = NaN;
        pH_L_D_20(end+1,1) = NaN;
        pH_L_S_20(end+1,1) = NaN;
        pH_L_D_corr_20(end+1,1) = NaN; 
        pH_L_S_corr_20(end+1,1) = NaN;
        pH_adj_20(end+1,1) = NaN;
        S_20(end+1,1) = NaN;
        pres_20(end+1,1) = NaN;
        temp_20(end+1,1) = NaN;
        S55_comp_float_data.Alk_pres(1:size(oxy(msk,i),1),end+1) = NaN;
    end   
end
% there are trailing 0s in the alkalinity results because the vectors are
% not equally long
Alk_LIAR(Alk_LIAR==0)=NaN;

S55_comp_float_data.Alk_LIAR = Alk_LIAR;


%Si ave from RAS SOFS8 2.8 umol/kg
%PO4 ave from RAS SOFS8 0.9
%NH4 set to 2 umol/kg
%H2S set to 0

% calculate fCO2 for float profiles with CO2SYS
%SOTS_float_data=[];
% I will use pCO2 for the floats, since I have pCO2 from Cape Grim
for i=1:length(Alk_LIAR_20)
    [DATA_L_D, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_L_D(i) = DATA_L_D(22); % uatm
    
    [DATA_L_D_20, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_20(i),1,3,S_20(i),temp_20(i),temp_20(i),20,20,2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_L_D_20m(i) = DATA_L_D_20(22); % uatm
    
    [DATA_L_S, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_S_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_L_S(i) = DATA_L_S(22); % uatm
   
    [DATA_L_D_corr, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_D_corr_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_L_D_corr(i) = DATA_L_D_corr(22); % uatm
    [DATA_L_S_corr, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_L_S_corr_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_L_S_corr(i) = DATA_L_S_corr(22); % uatm
  
    [DATA_pH_adj, HEADERS, NICEHEADERS]  = CO2SYS(Alk_LIAR_20(i),pH_adj_20(i),1,3,S_20(i),temp_20(i),temp_20(i),pres_20(i),pres_20(i),2.8,0.9,2,0,1,10,1,2,2);
    S55_comp_float_data.pCO2_pH_adj(i) = DATA_pH_adj(22); % uatm
  
    
    S55_comp_float_data.lat(i) = lat(i);
    S55_comp_float_data.lon(i) = lon(i);
    S55_comp_float_data.time(i) = time_float(i);
    S55_comp_float_data.TEMP_20(i) = temp_20(i);
    S55_comp_float_data.PSAL_20(i) = S_20(i);
    S55_comp_float_data.pH_L_D_20(i) = pH_L_D_20(i);
    S55_comp_float_data.pH_L_S_20(i) = pH_L_S_20(i);
    S55_comp_float_data.pH_L_D_20_corr(i) = pH_L_S_corr_20(i);
    S55_comp_float_data.pH_L_S_20_corr(i) = pH_L_S_corr_20(i);
    S55_comp_float_data.pH_adj_20(i) = pH_adj_20(i);
       
    S55_comp_float_data.Alk_LIAR_20(i) = Alk_LIAR_20(i);
    S55_comp_float_data.pres_20(i) = pres_20(i);
end


% read the ERA5 data for the area of the float
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript';
file2016 = [file_path '\2016_adaptor.mars.internal-1642655543.3866887-7060-10-8984e189-368a-4640-bd31-c0d4309a31ab.nc'];
file2017_18 = [file_path '\2017_18_adaptor.mars.internal-1642657799.05419-18345-16-60f5e400-895d-41d0-8cb9-66323e98beda.nc'];

Era2016=[];
Era2016.u_10 = ncread(file2016,'u10');
Era2016.v_10 = ncread(file2016, 'v10');
Era2016.t_2 = ncread(file2016, 't2m')-273.15; % Kelvin converted to Celsius
Era2016.msp = ncread(file2016, 'msl');
Era2016.EraTime = hours(ncread(file2016, 'time'))+ datetime(1900,1,1);
Era2016.EraLat = ncread(file2016, 'latitude');
Era2016.EraLon = ncread(file2016, 'longitude');

Era2017_18=[];
Era2017_18.u_10 = ncread(file2017_18,'u10');
Era2017_18.v_10 = ncread(file2017_18, 'v10');
Era2017_18.t_2 = ncread(file2017_18, 't2m')-273.15; % Kelvin converted to Celsius
Era2017_18.msp = ncread(file2017_18, 'msl');
Era2017_18.EraTime = hours(ncread(file2017_18, 'time'))+ datetime(1900,1,1);
Era2017_18.EraLat = ncread(file2017_18, 'latitude');
Era2017_18.EraLon = ncread(file2017_18, 'longitude');


% if float_ID == '5906623'
%     load('NCEP_SOTS_float.mat')
% elseif float_ID == '5906624'
%     load ('NCEP_S55_float.mat')
% end
% 
for i = 1:fs(2)
%for i = 156   
    if i<=31
        % find the ERA file index that is closest in time and space to the float
        data = Era2016;
        idx_T = knnsearch(datenum(data.EraTime(:)),datenum(time_float(i)));
        delta_time = (data.EraTime(idx_T)-time_float(i))/24;
        idx_lat = knnsearch(data.EraLat(:),lat(i));
        idx_lon = knnsearch(data.EraLon(:),lon(i));
    elseif i>31
        data=Era2017_18;
        % find the ERA file index that is closest in time and space to the float
        idx_T = knnsearch(datenum(data.EraTime(:)),datenum(time_float(i)));
        delta_time = (data.EraTime(idx_T)-time_float(i))/24;
        idx_lat = knnsearch(data.EraLat(:),lat(i));
        idx_lon = knnsearch(data.EraLon(:),lon(i));
    end
    
    if i<=31
        data = Era2016;
    elseif i>31
        data=Era2017_18;
    end

    wsp_f = sqrt(data.u_10(idx_lon,idx_lat,idx_T)^2 + data.v_10(idx_lon,idx_lat,idx_T)^2);
    %m s-2
    t_2_f = data.t_2(idx_lon,idx_lat,idx_T); %C
    msp_f = data.msp(idx_lon,idx_lat,idx_T); %Pa
    if abs(delta_time)>1
        wsp_f = NaN;
        t_2_f = NaN;
        msp_f = NaN;
    end
%         
%     if isnan(wsp_f)
%         wsp_f = float_ncep.wnd_ncep_f(i);
%     end
%     
%     if isnan(t_2_f)
%         t_2_f = float_ncep.airT_ncep_f(i);
%     end
%         
%     if isnan(msp_f)
%         msp_f = float_ncep.mslp_ncep_f(i);
%     end
        
    
    % water vapour pressure 6.11 * 10^((7.5*(T))/(237.3+(T)))
    % Temp in degrees C 
    pH20_f = 6.11 * 10^((7.5*t_2_f)/(273.3+t_2_f)); %mbar = hPa

    % convert xCO2 in dry air to pCO2
    idx_C = knnsearch(datenum(CGdata_interp.spl_time(:)),datenum(time_float(i)));
    %uatm
    pCO2_uatm = CGdata_interp.ppm_spl(idx_C) * ((msp_f/100 - pH20_f))*(9.8692326671601*10^-4); %uatm
    
    % now combine this with atmospheric fCO2
    DpCO2_L_D = S55_comp_float_data.pCO2_L_D(i)-pCO2_uatm;
    DpCO2_L_D_20m = S55_comp_float_data.pCO2_L_D_20m(i)-pCO2_uatm;
    DpCO2_L_S = S55_comp_float_data.pCO2_L_S(i)-pCO2_uatm;
    DpCO2_L_D_corr = S55_comp_float_data.pCO2_L_D_corr(i)-pCO2_uatm;
    DpCO2_L_S_corr = S55_comp_float_data.pCO2_L_S_corr(i)-pCO2_uatm;
    DpCO2_pH_adj = S55_comp_float_data.pCO2_pH_adj(i)-pCO2_uatm;

    
    [F_CO2_float_L_D]=FCO2_CWE(DpCO2_L_D,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_D_20m]=FCO2_CWE(DpCO2_L_D_20m,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_S]=FCO2_CWE(DpCO2_L_S,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_D_corr]=FCO2_CWE(DpCO2_L_D_corr,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_L_S_corr]=FCO2_CWE(DpCO2_L_S_corr,temp_20(i),S_20(i),wsp_f);
    [F_CO2_float_pH_adj]=FCO2_CWE(DpCO2_pH_adj,temp_20(i),S_20(i),wsp_f);
   
   
    S55_comp_float_data.flux_L_D(i) = F_CO2_float_L_D;
    S55_comp_float_data.flux_L_D_20m(i) = F_CO2_float_L_D_20m;
    S55_comp_float_data.flux_L_S(i) = F_CO2_float_L_S;
    S55_comp_float_data.flux_L_D_corr(i) = F_CO2_float_L_D_corr;
    S55_comp_float_data.flux_L_S_corr(i) = F_CO2_float_L_S_corr;
    S55_comp_float_data.flux_pH_adj(i) = F_CO2_float_pH_adj;
    
   
    S55_comp_float_data.pCO2_uatm(i) = pCO2_uatm;
    S55_comp_float_data.CG_xCO2_uatm(i) = CGdata_interp.ppm_spl(idx_C);
    S55_comp_float_data.wsp(i) = wsp_f;
    S55_comp_float_data.t_2(i) = t_2_f;
    S55_comp_float_data.msp(i) = msp_f;
end


%%%% now I just want 2016
S55_comp_float_data.Month = month(S55_comp_float_data.time);
S55_comp_float_data.year = year(S55_comp_float_data.time);

S55_comp_float_data.data_comp_2016 = table(S55_comp_float_data.Month(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.time(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.pH_L_D_20(S55_comp_float_data.year==2016)',S55_comp_float_data.pH_L_S_20(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.pH_L_D_20_corr(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.pH_L_S_20_corr(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.pH_adj_20(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.Alk_LIAR_20(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.flux_L_D(S55_comp_float_data.year==2016)', S55_comp_float_data.flux_L_S(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.flux_L_D_corr(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.flux_L_S_corr(S55_comp_float_data.year==2016)',...
        S55_comp_float_data.flux_pH_adj(S55_comp_float_data.year==2016)',...    
         'VariableNames',{'Month','Time','pH_LD_20','pH_LS_20',...
        'pH_LD_20_corr','pH_LS_20_corr','pH_adj_20','Alk_LIAR_20',...
        'flux_LD','flux_LS','flux_LD_corr','flux_LS_corr','flux_pH_adj'});

S55_comp_float_data.y2016.pH_LD_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.pH_LD_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.pH_LS_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.pH_LS_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.pH_LD_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.pH_LD_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.pH_LS_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.pH_LS_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.pH_adj_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.pH_adj_20,[],@(x)mean(x,'omitnan'));


S55_comp_float_data.y2016.Alk_LIAR_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.Alk_LIAR_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.flux_LD_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.flux_LD,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.flux_LS_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.flux_LS,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.flux_LD_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.flux_LD_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.flux_LS_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.flux_LS_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2016.flux_pH_adj_mo_ave = accumarray(S55_comp_float_data.data_comp_2016.Month,S55_comp_float_data.data_comp_2016.flux_pH_adj,[],@(x)mean(x,'omitnan'));

S55_comp_float_data.y2016.mo_ave_month = unique(S55_comp_float_data.data_comp_2016.Month);


%%%% now I just want 2017

S55_comp_float_data.data_comp_2017 = table(S55_comp_float_data.Month(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.time(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.pH_L_D_20(S55_comp_float_data.year==2017)',S55_comp_float_data.pH_L_S_20(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.pH_L_D_20_corr(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.pH_L_S_20_corr(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.pH_adj_20(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.Alk_LIAR_20(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.flux_L_D(S55_comp_float_data.year==2017)', S55_comp_float_data.flux_L_S(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.flux_L_D_corr(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.flux_L_S_corr(S55_comp_float_data.year==2017)',...
        S55_comp_float_data.flux_pH_adj(S55_comp_float_data.year==2017)',...    
         'VariableNames',{'Month','Time','pH_LD_20','pH_LS_20',...
        'pH_LD_20_corr','pH_LS_20_corr','pH_adj_20','Alk_LIAR_20',...
        'flux_LD','flux_LS','flux_LD_corr','flux_LS_corr','flux_pH_adj'});

S55_comp_float_data.y2017.pH_LD_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.pH_LD_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.pH_LS_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.pH_LS_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.pH_LD_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.pH_LD_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.pH_LS_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.pH_LS_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.pH_adj_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.pH_adj_20,[],@(x)mean(x,'omitnan'));


S55_comp_float_data.y2017.Alk_LIAR_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.Alk_LIAR_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.flux_LD_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.flux_LD,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.flux_LS_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.flux_LS,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.flux_LD_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.flux_LD_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.flux_LS_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.flux_LS_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2017.flux_pH_adj_mo_ave = accumarray(S55_comp_float_data.data_comp_2017.Month,S55_comp_float_data.data_comp_2017.flux_pH_adj,[],@(x)mean(x,'omitnan'));

S55_comp_float_data.y2017.mo_ave_month = unique(S55_comp_float_data.data_comp_2017.Month);

%%%% now I just want 2018

S55_comp_float_data.data_comp_2018 = table(S55_comp_float_data.Month(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.time(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.pH_L_D_20(S55_comp_float_data.year==2018)',S55_comp_float_data.pH_L_S_20(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.pH_L_D_20_corr(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.pH_L_S_20_corr(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.pH_adj_20(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.Alk_LIAR_20(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.flux_L_D(S55_comp_float_data.year==2018)', S55_comp_float_data.flux_L_S(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.flux_L_D_corr(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.flux_L_S_corr(S55_comp_float_data.year==2018)',...
        S55_comp_float_data.flux_pH_adj(S55_comp_float_data.year==2018)',...    
         'VariableNames',{'Month','Time','pH_LD_20','pH_LS_20',...
        'pH_LD_20_corr','pH_LS_20_corr','pH_adj_20','Alk_LIAR_20',...
        'flux_LD','flux_LS','flux_LD_corr','flux_LS_corr','flux_pH_adj'});

S55_comp_float_data.y2018.pH_LD_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.pH_LD_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.pH_LS_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.pH_LS_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.pH_LD_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.pH_LD_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.pH_LS_20_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.pH_LS_20_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.pH_adj_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.pH_adj_20,[],@(x)mean(x,'omitnan'));


S55_comp_float_data.y2018.Alk_LIAR_20_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.Alk_LIAR_20,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.flux_LD_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.flux_LD,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.flux_LS_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.flux_LS,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.flux_LD_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.flux_LD_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.flux_LS_corr_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.flux_LS_corr,[],@(x)mean(x,'omitnan'));
S55_comp_float_data.y2018.flux_pH_adj_mo_ave = accumarray(S55_comp_float_data.data_comp_2018.Month,S55_comp_float_data.data_comp_2018.flux_pH_adj,[],@(x)mean(x,'omitnan'));

S55_comp_float_data.y2018.mo_ave_month = unique(S55_comp_float_data.data_comp_2018.Month);



clearvars -except S55_comp_float_data pH_ref CGdata_interp float_ID

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('S55_comp_float_data.mat')