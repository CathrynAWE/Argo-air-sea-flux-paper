%%%% trying to work out why float raw data pH correction with LIPHR deep is
%%%% a better fit than LIPHR shallow, 55S float


load('S55_float_data.mat')

addpath(genpath('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0'))


figure()
plot(S55_float_data.time,pH_ref.pH_D,'ok')
hold on
plot(S55_float_data.time,pH_ref.pH_S,'*r')
xlabel('profile date')
ylabel('pH at ref depth')
legend('1500m','950m')


figure()
plot(S55_float_data.time,S55_float_data.pH_L_D_20_corr,'ok')
hold on
plot(S55_float_data.time,S55_float_data.pH_L_S_20_corr,'*r')
xlabel('profile date')
ylabel('average pH of top 20m')
legend('LIPHR deep','LIPHR shallow')


float_ID = '5906624' % 55S float
% %%%%% change structure name!
% S55_float_data = [];
% 
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

%%% calculate pH based on float S, O2 and temp
for i = 1:fs(2)
%for i =1:3
    % mask out the NaNs in the float data, for the Alklinity script
%     msk = ~isnan(oxy(:,i));
    lon_c = zeros(size(oxy(:,i),1),1);
    lon_c(:) = lon(i);
    lat_c = zeros(size(oxy(:,i),1),1);
    lat_c(:) = lat(1);
    Coor = [lon_c lat_c pres(:,i)];

    Meas = [S(:,i) oxy(:,i) temp(:,i)];
    
    % define the measurement IDs (salinity, oxygen, temperature)
    MeasID = [1 6 7];

    % run the LIPHR code to estimate pH
    disp(i)
    pH_corr.pH_LIPHR(:,i) = LIPHR(Coor, Meas, MeasID,'OAAdjustTF',false);
    pH_corr.pH_LIPHR_OA(:,i) = LIPHR(Coor, Meas, MeasID, 'EstDates', 2021.5,'OAAdjustTF',true);
                  
end

%%% find the pH at the deep reference depth (1500m)
%%%%% set the reference depth (1500m as per Williams, 2017)
%%%%% I will use this for the LD and WD corrected profiles
ref_depth_D = 1500;
%%% set the tolerance to the reference depth
ref_depth_t_D = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i= 1:size(pH,2)
    pres_ref_id_D = pres(:,i)>=(ref_depth_D-ref_depth_t_D) & pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
    pres_D(i) = mean(pres(pres_ref_id_D,i),'omitnan');
    %%% find the pH at the reference depths
    pH_corr.pH_D(i) = mean(pH(pres_ref_id_D,i),'omitnan'); % raw pH at 1500m
    pH_corr.pH_LIPHR_D(i)=mean(pH_corr.pH_LIPHR(pres_ref_id_D,i),'omitnan'); % LIPHR estimate at 1500m
    pH_corr.pH_LIPHR_D_OA(i)=mean(pH_corr.pH_LIPHR_OA(pres_ref_id_D,i),'omitnan'); % LIPHR estimate at 1500m, with OA adjustment
    pH_corr.pH_D_corr(i) = mean(S55_float_data.pH_LIR_Deep(pres_ref_id_D,i),'omitnan'); % pH corrected at 1500m
end


%%% find the pH at the deep reference depth (950m)
%%%%% I will use this for the LS profiles
ref_depth_S = 950;
%%% set the tolerance to the reference depth
ref_depth_t_S = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i= 1:size(pH,2)
    pres_ref_id_S = pres(:,i)>=(ref_depth_S-ref_depth_t_S) & pres(:,i)<=(ref_depth_S+ref_depth_t_S); 
    pres_S(i) = mean(pres(pres_ref_id_S,i),'omitnan');
    %%% find the pH at the reference depths
    pH_corr.pH_S(i) = mean(pH(pres_ref_id_S,i),'omitnan') ;% raw pH at 950m
    pH_corr.pH_LIPHR_S(i)=mean(pH_corr.pH_LIPHR(pres_ref_id_S,i),'omitnan'); % LIPHR estimate at 950m
    pH_corr.pH_LIPHR_S_OA(i)=mean(pH_corr.pH_LIPHR_OA(pres_ref_id_S,i),'omitnan'); % LIPHR estimate at 950m, with OA adjustment
    pH_corr.pH_S_corr(i) = mean(S55_float_data.pH_LIR_Shallow(pres_ref_id_S,i),'omitnan'); % pH corrected at 950m
end


%%% for comparison, I now want to find the ph at both depths, once in the
%%% LD corrected profiles and once in the LS corrected profiles so I can
%%% plot one vs the other, like I did for raw pH and LIPHR estimates
for i= 1:size(pH,2)
    pres_ref_id_D = pres(:,i)>=(ref_depth_D-ref_depth_t_D) & pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
    pres_ref_id_S = pres(:,i)>=(ref_depth_S-ref_depth_t_S) & pres(:,i)<=(ref_depth_D+ref_depth_t_S); 
    %%% find the pH at the reference depths
    pH_corr.pH_D_corr_LS(i) = mean(S55_float_data.pH_LIR_Shallow(pres_ref_id_D,i),'omitnan'); %this is the pH at 1500m in the LS corrected profiles
    pH_corr.pH_S_corr_LD(i) = mean(S55_float_data.pH_LIR_Deep(pres_ref_id_S,i),'omitnan'); %this is the pH at 950m in the LD corrected profiles
end

cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper')
save('pH-corr-cal_S55.mat','pH_corr')

% load('pH-corr-cal.mat')

figure()
plot(pH_corr.pH_LIPHR_D,pH_corr.pH_D,'.k')
hold on
plot(pH_corr.pH_LIPHR_D_OA,pH_corr.pH_D,'.g')
xlabel('LIPHR pH at 1500m')
ylabel('float raw pH at 1500m')
legend('without OA', 'with OA')

lr_LIHPR_D = fitlm(pH_corr.pH_LIPHR_D,pH_corr.pH_D)
lr_LIPHR_D_OA = fitlm(pH_corr.pH_LIPHR_D_OA,pH_corr.pH_D)


figure()
plot(pH_corr.pH_LIPHR_S,pH_corr.pH_S,'.r')
hold on
plot(pH_corr.pH_LIPHR_S_OA,pH_corr.pH_S,'.g')
xlabel('LIPHR pH at 950m')
ylabel('float raw pH at 950m')
legend('without OA', 'with OA')


lr_LIHPR_S = fitlm(pH_corr.pH_LIPHR_S,pH_corr.pH_S)
lr_LIPHR_S_OA = fitlm(pH_corr.pH_LIPHR_S_OA,pH_corr.pH_S)



figure()
plot(pH_corr.pH_D,pH_corr.pH_S,'.','MarkerSize',6)
xlabel('pH at 1500m')
ylabel('pH at 950m')
hold on
plot(pH_corr.pH_LIPHR_D,pH_corr.pH_LIPHR_S,'.r','MarkerSize',12)
plot(pH_corr.pH_LIPHR_D_OA,pH_corr.pH_LIPHR_S_OA,'.k','MarkerSize',12)
plot(pH_corr.pH_D_corr,pH_corr.pH_S_corr_LD,'^k','MarkerSize',4)
plot(pH_corr.pH_D_corr_LS,pH_corr.pH_S_corr,'^g','MarkerSize',4)
legend('raw','LIPHR no OA','LIPHR with OA','corrected LD','corrected LS')

mdl_raw=fitlm(pH_corr.pH_D,pH_corr.pH_S);
mdl_LIPHR=fitlm(pH_corr.pH_LIPHR_D,pH_corr.pH_LIPHR_S);
mdl_LIPHR_OA=fitlm(pH_corr.pH_LIPHR_D_OA,pH_corr.pH_LIPHR_S_OA);
mdl_pH_corr_LD=fitlm(pH_corr.pH_D_corr,pH_corr.pH_S_corr_LD);
mdl_pH_corr_LS=fitlm(pH_corr.pH_D_corr_LS,pH_corr.pH_S_corr);

figure()
plot(mdl_raw)
hold on
plot(mdl_LIPHR)
plot(mdl_LIPHR_OA)
plot(mdl_pH_corr_LD)
plot(mdl_pH_corr_LS)

plot(pH_corr.pH_D,pH_corr.pH_S,'.','MarkerSize',6)
plot(pH_corr.pH_LIPHR_D,pH_corr.pH_LIPHR_S,'.r','MarkerSize',12)
plot(pH_corr.pH_LIPHR_D_OA,pH_corr.pH_LIPHR_S_OA,'.k','MarkerSize',12)
plot(pH_corr.pH_D_corr,pH_corr.pH_S_corr_LD,'^k','MarkerSize',4)
plot(pH_corr.pH_D_corr_LS,pH_corr.pH_S_corr,'^g','MarkerSize',4)
xlabel('pH at 1500m')
ylabel('pH at 950m')
legend('raw','LIPHR','corrected LD','corrected LS')


offset=S55_float_data.pH_L_D_20_corr-S55_float_data.pH_L_S_20_corr;

figure()
plot(S55_float_data.time,offset,'.')



