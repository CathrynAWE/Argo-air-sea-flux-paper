%%%% trying to work out why float raw data pH correction with LIPHR deep is
%%%% a better fit than LIPHR shallow


load('SOTS_float_data_TEMP.mat')

addpath(genpath('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper\seawater_ver3_0'))
% 
% figure()
% plot(SOTS_float_data_TEMP.pH(SOTS_float_data_TEMP.pres<=20),SOTS_float_data_TEMP.pH_LIR_Deep(SOTS_float_data_TEMP.pres<=20),'.k')
% xlabel('raw pH top 20m')
% ylabel('pH LIPHR Deep top 20m')
% ylim([8.01 8.15])
% 
% figure()
% plot(SOTS_float_data_TEMP.pH(SOTS_float_data_TEMP.pres<=20),SOTS_float_data_TEMP.pH_LIR_Shallow(SOTS_float_data_TEMP.pres<=20),'.r')
% xlabel('raw pH top 20m')
% ylabel('pH LIPHR shallow top 20m')
% ylim([8.01 8.15])


%%% float pH vs SOLACE CTD cast pH
%%% with TS diagram
load('CTD_data.mat')
% what if I limit the data points to compare to those where CTD data
% and float data are within 1 or 2 days and with temp and psal limits that
% equate to those set for well mixed waters (mixed layer depth
% definitions), to make sure we only compare data from the same water
% masses. Temp threshold = 0.3C (temp QC report), psal threshold = 0.03psu
% (salinity QC report)

% time threshold in days
d = 2;
% temp threshold in C
t = 0.3;
% psal threshold in psu
ps = 0.03;
% pres threshold in dbar
p = 5;

% create the empty table to hold the final subset, with vartypes and names
% of the original CTD table
vart = varfun(@class,CTD_data.raw_data,'OutputFormat','cell');
varname = CTD_data.raw_data.Properties.VariableNames;
s=[1 size(CTD_data.raw_data,2)];

CTD_pH_subset=table('Size',s,'VariableNames', varname,'VariableTypes',vart);


for i = 1:size(SOTS_float_data_TEMP.pH,2)
    %subset by time threshold
    CTD_temporal_sub = CTD_data.raw_data(abs(datenum(CTD_data.raw_data.date)-datenum(SOTS_float_data_TEMP.time(i)))<=d,:);
    if ~isempty(CTD_temporal_sub)
        pH_idx(i)=1;
        for j = 1:size(SOTS_float_data_TEMP.pH,1)
            % find the CTD data the fit within the pressure threshold of each
            % float measurement
            pres_sub = CTD_temporal_sub(abs(CTD_temporal_sub.Depth - SOTS_float_data_TEMP.pres(j,i))<=p,:);
            if ~isempty(pres_sub)
                pres_idx(j)=1;
                %subset by temp threshold
                temp_sub = pres_sub(abs(pres_sub.T_insitu-SOTS_float_data_TEMP.TEMP(j,i))<=t,:);
                if ~isempty(temp_sub)
                    temp_idx(j)=1;
                    %subset by psal threshold
                    psal_sub = temp_sub(abs(temp_sub.Salinity_CTD-SOTS_float_data_TEMP.psal(j,i))<=ps,:);
                    if ~isempty(psal_sub)
                        psal_idx(j)=1;
                        %now add this to the CTD_pH_subset to be plotted
                        %against float data
                        for h=1:height(psal_sub)
                            CTD_pH_subset(end+1,:) = psal_sub(h,:);
                        end
                    end
                end
            end
        end
        if sum(pres_idx)+sum(temp_idx)+sum(psal_idx)==0
            pH_idx(i)=0;
        end
    else
        pH_idx(i)=0;
    end
end

% clear out the first empty line
CTD_pH_subset(1,:)=[];

cols = lines; 

figure()
subplot(1,2,1)
plot(SOTS_float_data_TEMP.pH_LD_corr(:,pH_idx==1), SOTS_float_data_TEMP.pres(:,pH_idx==1), '.', 'MarkerSize', 10,'color',cols(3,:))
hold on
plot(SOTS_float_data_TEMP.pH_LS_corr(:,pH_idx==1), SOTS_float_data_TEMP.pres(:,pH_idx==1), '.g', 'MarkerSize', 10,'color',cols(5,:))

plot(CTD_pH_subset.pH, CTD_pH_subset.Depth,'.b', 'MarkerSize', 20)
hold off
xlabel('pH total scale')
ylabel('Depth dbar')
legend('\color{orange}float LD corr','\color{green}float LS corr','\color{blue}CTD','Orientation','horizontal','Location','bestoutside')
set(gca, 'YDir','reverse')
ylim([0 1600])

subplot(1,2,2)
plot(SOTS_float_data_TEMP.TEMP(:,pH_idx==1), SOTS_float_data_TEMP.psal(:,pH_idx==1), '^k', 'MarkerSize', 4)
hold on
plot(CTD_pH_subset.T_insitu, CTD_pH_subset.Salinity_CTD,'.b', 'MarkerSize', 20)
hold off
xlabel('Temp')
ylabel('PSAL')
legend('\color{black}float','\color{blue}CTD','Orientation','horizontal','Location','bestoutside')



figure()
plot(SOTS_float_data_TEMP.time,pH_ref.pH_D,'ok')
hold on
plot(SOTS_float_data_TEMP.time,pH_ref.pH_S,'*r')
xlabel('profile date')
ylabel('pH at ref depth')
legend('1500m','950m')


figure()
plot(SOTS_float_data_TEMP.time,SOTS_float_data_TEMP.pH_L_D_20_corr,'ok')
hold on
plot(SOTS_float_data_TEMP.time,SOTS_float_data_TEMP.pH_L_S_20_corr,'*r')
xlabel('profile date')
ylabel('average pH of top 20m')
legend('LIPHR deep','LIPHR shallow')


float_ID = '5906623' % SOTS float
% float_ID = '5906624' % 55S float
%%%%% change structure name!
% SOTS_float_data_TEMP = [];

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
% 
% %%%% clear out the bad profiles and all their successive calculations
% if float_ID == '5906623'
%     pH(:,10:11) = NaN;
%     pH(:,23) = NaN;
% elseif float_ID == '5906624'
%     pH(:,1:3) = NaN;
% end
% 
% %%% calculate pH based on float S, O2 and temp
% for i = 1:fs(2)
% %for i =1:3
%     % mask out the NaNs in the float data, for the Alklinity script
% %     msk = ~isnan(oxy(:,i));
%     lon_c = zeros(size(oxy(:,i),1),1);
%     lon_c(:) = lon(i);
%     lat_c = zeros(size(oxy(:,i),1),1);
%     lat_c(:) = lat(1);
%     Coor = [lon_c lat_c pres(:,i)];
% 
%     Meas = [S(:,i) oxy(:,i) temp(:,i)];
%     
%     % define the measurement IDs (salinity, oxygen, temperature)
%     MeasID = [1 6 7];
% 
%     % run the LIPHR code to estimate pH
%     disp(i)
%     pH_corr.pH_LIPHR(:,i) = LIPHR(Coor, Meas, MeasID,'OAAdjustTF',false);
%     pH_corr.pH_LIPHR_OA(:,i) = LIPHR(Coor, Meas, MeasID, 'EstDates', 2021.5,'OAAdjustTF',true);
%                   
% end
% 
% %%% find the pH at the deep reference depth (1500m)
% %%%%% set the reference depth (1500m as per Williams, 2017)
% %%%%% I will use this for the LD and WD corrected profiles
% ref_depth_D = 1500;
% %%% set the tolerance to the reference depth
% ref_depth_t_D = 50;
% %%% find the index of each float profile where the depth is closest to the
% %%% reference depth
% for i= 1:size(pH,2)
%     pres_ref_id_D = pres(:,i)>=(ref_depth_D-ref_depth_t_D) & pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
%     pres_D(i) = mean(pres(pres_ref_id_D,i),'omitnan');
%     %%% find the pH at the reference depths
%     pH_corr.pH_D(i) = mean(pH(pres_ref_id_D,i),'omitnan'); % raw pH at 1500m
%     pH_corr.pH_LIPHR_D(i)=mean(pH_corr.pH_LIPHR(pres_ref_id_D,i),'omitnan'); % LIPHR estimate at 1500m
%     pH_corr.pH_LIPHR_D_OA(i)=mean(pH_corr.pH_LIPHR_OA(pres_ref_id_D,i),'omitnan'); % LIPHR estimate at 1500m, with OA adjustment
%     pH_corr.pH_D_corr(i) = mean(SOTS_float_data_TEMP.pH_LIR_Deep(pres_ref_id_D,i),'omitnan'); % pH corrected at 1500m
% end
% 
% 
% %%% find the pH at the deep reference depth (950m)
% %%%%% I will use this for the LS profiles
% ref_depth_S = 950;
% %%% set the tolerance to the reference depth
% ref_depth_t_S = 50;
% %%% find the index of each float profile where the depth is closest to the
% %%% reference depth
% for i= 1:size(pH,2)
%     pres_ref_id_S = pres(:,i)>=(ref_depth_S-ref_depth_t_S) & pres(:,i)<=(ref_depth_S+ref_depth_t_S); 
%     pres_S(i) = mean(pres(pres_ref_id_S,i),'omitnan');
%     %%% find the pH at the reference depths
%     pH_corr.pH_S(i) = mean(pH(pres_ref_id_S,i),'omitnan') ;% raw pH at 950m
%     pH_corr.pH_LIPHR_S(i)=mean(pH_corr.pH_LIPHR(pres_ref_id_S,i),'omitnan'); % LIPHR estimate at 950m
%     pH_corr.pH_LIPHR_S_OA(i)=mean(pH_corr.pH_LIPHR_OA(pres_ref_id_S,i),'omitnan'); % LIPHR estimate at 950m, with OA adjustment
%     pH_corr.pH_S_corr(i) = mean(SOTS_float_data_TEMP.pH_LIR_Shallow(pres_ref_id_S,i),'omitnan'); % pH corrected at 950m
% end
% 
% 
% %%% for comparison, I now want to find the ph at both depths, once in the
% %%% LD corrected profiles and once in the LS corrected profiles so I can
% %%% plot one vs the other, like I did for raw pH and LIPHR estimates
% for i= 1:size(pH,2)
%     pres_ref_id_D = pres(:,i)>=(ref_depth_D-ref_depth_t_D) & pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
%     pres_ref_id_S = pres(:,i)>=(ref_depth_S-ref_depth_t_S) & pres(:,i)<=(ref_depth_D+ref_depth_t_S); 
%     %%% find the pH at the reference depths
%     pH_corr.pH_D_corr_LS(i) = mean(SOTS_float_data_TEMP.pH_LIR_Shallow(pres_ref_id_D,i),'omitnan'); %this is the pH at 1500m in the LS corrected profiles
%     pH_corr.pH_S_corr_LD(i) = mean(SOTS_float_data_TEMP.pH_LIR_Deep(pres_ref_id_S,i),'omitnan'); %this is the pH at 950m in the LD corrected profiles
% end

% cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper')
% save('pH-corr-cal.mat','pH_corr')

load('pH-corr-cal.mat')

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
mdl_pH_corr_LD=fitlm(pH_corr.pH_D_corr,pH_corr.pH_S_corr_LD);
mdl_pH_corr_LS=fitlm(pH_corr.pH_D_corr_LS,pH_corr.pH_S_corr);

figure()
plot(mdl_raw)
hold on
plot(mdl_LIPHR)
plot(mdl_pH_corr_LD)
plot(mdl_pH_corr_LS)
plot(pH_corr.pH_D,pH_corr.pH_S,'.','MarkerSize',6)
plot(pH_corr.pH_LIPHR_D,pH_corr.pH_LIPHR_S,'.r','MarkerSize',12)
plot(pH_corr.pH_D_corr,pH_corr.pH_S_corr_LD,'^k','MarkerSize',4)
plot(pH_corr.pH_D_corr_LS,pH_corr.pH_S_corr,'^g','MarkerSize',4)
xlabel('pH at 1500m')
ylabel('pH at 950m')
legend('raw','LIPHR','corrected LD','corrected LS')


offset=SOTS_float_data_TEMP.pH_L_D_20_corr-SOTS_float_data_TEMP.pH_L_S_20_corr;

figure()
plot(SOTS_float_data_TEMP.time,offset,'.')




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% all float pH vs CTD casts at SOTS
load('CTD_data.mat')

figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data_TEMP.pH_LIR_Deep(:,month(SOTS_float_data_TEMP.time)==i),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==i),'^r', 'MarkerSize',2)
    hold on
    plot(SOTS_float_data_TEMP.pH_LIR_Shallow(:,month(SOTS_float_data_TEMP.time)==i),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==i),'^g', 'MarkerSize',2)
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==i), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==i),'.b','MarkerSize',10)
    plot(pH_corr.pH_LIPHR(:,month(SOTS_float_data_TEMP.time)==i),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==i),'*m', 'MarkerSize',2)
    if i ==1
        xlabel('pH total scale - Jan')
    elseif i ==2
        xlabel('pH total scale - Feb')
    elseif i==3
        xlabel('pH total scale - Mar')
    elseif i==4
        xlabel('pH total scale - Apr')
    elseif i==5
        xlabel('pH total scale - May')
    elseif i==6
        xlabel('pH total scale - June')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
    
end
legend('\color{red}float LD','\color{green}float LS','\color{blue}CTD','\color{magenta}float LIPHR','Location','best')

figure()
for i = 1:6
    subplot(3,2,i)
    plot(SOTS_float_data_TEMP.pH_LD_corr(:,month(SOTS_float_data_TEMP.time)==(i+6)),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==(i+6)),'^r', 'MarkerSize',2)% akin to pH estimated from DIC and Alk (as per Takeshita)
    hold on
    plot(SOTS_float_data_TEMP.pH_LS_corr(:,month(SOTS_float_data_TEMP.time)==(i+6)),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==(i+6)),'^g', 'MarkerSize',2)% akin to pH estimated from DIC and Alk (as per Takeshita)
%     plot(SOTS_float_data.pH_Williams_Deep(:,month(SOTS_float_data.time)==(i+6)),SOTS_float_data.pres(:,month(SOTS_float_data.time)==(i+6)),'^c', 'MarkerSize',2)
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==(i+6)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==(i+6)),'.b','MarkerSize',10)
    plot(pH_corr.pH_LIPHR(:,month(SOTS_float_data_TEMP.time)==i+6),SOTS_float_data_TEMP.pres(:,month(SOTS_float_data_TEMP.time)==i+6),'*m', 'MarkerSize',2)
    if i ==1
        xlabel('pH total scale - July')
    elseif i ==2
        xlabel('pH total scale - Aug')
    elseif i==3
        xlabel('pH total scale - Sep')
    elseif i==4
        xlabel('pH total scale - Oct')
    elseif i==5
        xlabel('pH total scale - Nov')
    elseif i==6
        xlabel('pH total scale - Dec')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
    
end
legend('\color{red}float LD','\color{green}float LS','\color{blue}CTD','\color{magenta}float LIPHR','Location','best')


%%%%%%%%%%%%%
%%%% now I want to see what the LIPHR linear regressions produce with CTD
%%%% oxygen, salinity and temp (like I use for the floats)
%%%% the CTD data is available in different formats and the csv files are
%%%% not all the same, so I have to patch the data together

CTD_data_csv = readtable('CTD_data_ox.csv');
CTD_data_nc = readtable('CTD_nc_data_ox.csv');
w=width(CTD_data_csv);
l=height(CTD_data_csv)+height(CTD_data_nc);
CTD_data_ox=array2table(zeros(l,w));
names = CTD_data_csv.Properties.VariableNames;
CTD_data_ox.Properties.VariableNames = names;

CTD_data_ox.SURVEY_NAME=[CTD_data_csv.SURVEY_NAME;CTD_data_nc.SURVEY_NAME];
CTD_data_ox.STATION=[CTD_data_csv.STATION;CTD_data_nc.STATION];
CTD_data_ox.STATION=[CTD_data_csv.STATION;CTD_data_nc.STATION];
CTD_data_ox.START_TIME=[CTD_data_csv.START_TIME;CTD_data_nc.START_TIME];
CTD_data_ox.END_TIME=[CTD_data_csv.END_TIME;CTD_data_nc.END_TIME];
CTD_data_ox.MIN_DEPTH=[CTD_data_csv.MIN_DEPTH;CTD_data_nc.MIN_DEPTH];
CTD_data_ox.MAX_DEPTH=[CTD_data_csv.MAX_DEPTH;CTD_data_nc.MAX_DEPTH];
CTD_data_ox.BOTTOM_LAT=[CTD_data_csv.BOTTOM_LAT;CTD_data_nc.BOTTOM_LAT];
CTD_data_ox.BOTTOM_LON=[CTD_data_csv.BOTTOM_LON;CTD_data_nc.BOTTOM_LON];
CTD_data_ox.START_LAT=[CTD_data_csv.START_LAT;CTD_data_nc.START_LAT];
CTD_data_ox.START_LON=[CTD_data_csv.START_LON;CTD_data_nc.START_LON];
CTD_data_ox.PRESSURE=[CTD_data_csv.PRESSURE;CTD_data_nc.PRESSURE];
CTD_data_ox.OXYGEN=[CTD_data_csv.OXYGEN;CTD_data_nc.OXYGEN];
CTD_data_ox.OXYGEN_QC=[CTD_data_csv.OXYGEN_QC;CTD_data_nc.OXYGEN_QC];
CTD_data_ox.TEMPERATURE=[CTD_data_csv.TEMPERATURE;CTD_data_nc.TEMPERATURE];
CTD_data_ox.TEMPERATURE_QC=[CTD_data_csv.TEMPERATURE_QC;CTD_data_nc.TEMPERATURE_QC];
CTD_data_ox.SALINITY=[CTD_data_csv.SALINITY;CTD_data_nc.SALINITY];
CTD_data_ox.SALINITY_QC=[CTD_data_csv.SALINITY_QC;CTD_data_nc.SALINITY_QC];
CTD_data_ox.BOTTOM_LAT(isnan(CTD_data_ox.BOTTOM_LAT))=CTD_data_ox.START_LAT(isnan(CTD_data_ox.BOTTOM_LAT));

Coor_ox = [CTD_data_ox.BOTTOM_LON, CTD_data_ox.BOTTOM_LAT, CTD_data_ox.PRESSURE];
Meas_ox = [CTD_data_ox.SALINITY, CTD_data_ox.OXYGEN, CTD_data_ox.TEMPERATURE];% I am using the sensor umol/L data, but could use the umol/kg that I calculated
% define the measurement IDs (salinity, oxygen, temperature)
MeasID = [1 6 7];


CTD_data_ox.pH_LIPHR = LIPHR(Coor_ox, Meas_ox, MeasID,'OAAdjustTF',false,'PerKgSwTF',false);%oxygen sensor data in umol/L
CTD_data_ox.pH_LIPHR_DIC_Alk = LIPHR(Coor_ox, Meas_ox, MeasID,'OAAdjustTF',false,'PerKgSwTF',false,'pHCalcTF',true);%oxygen sensor data in umol/L, and turned on pHCalcTF to get a better alignment with pH calculated from DIC and alk

estdate = years(CTD_data_ox.START_TIME - datetime(2000,1,1)) + 2000;
CTD_data_ox.pH_LIPHR_OA = LIPHR(Coor_ox, Meas_ox, MeasID, 'EstDates', estdate,'OAAdjustTF', true,'PerKgSwTF',false); %oxygen sensor data in umol/L
CTD_data_ox.pH_LIPHR_OA_DIC_Alk = LIPHR(Coor_ox, Meas_ox, MeasID, 'EstDates', estdate,'OAAdjustTF', true,'PerKgSwTF',false,'pHCalcTF',true); %oxygen sensor data in umol/L, and turned on pHCalcTF to get a better alignment with pH calculated from DIC and alk


%%%%%%%%%%%%%%%%%%%%%%
%%%%% now I want to find the CTD sensor data for each of the bottle data
%%%% first for each line in the CTD_data, create the subset in the
%%%% CTD_data_ox set that holds the right voyage and CTD cast
CTD_data.raw_data_plusox = CTD_data.raw_data;


for i = 1:height(CTD_data.raw_data)
    dummy = CTD_data_ox(strcmp(lower(CTD_data_ox.SURVEY_NAME), lower(CTD_data.raw_data.Sample(i))) & CTD_data_ox.STATION == CTD_data.raw_data.CTD(i),:);
    if ~isempty(dummy)
        idx = knnsearch(dummy.PRESSURE, CTD_data.raw_data.Depth(i));
        if ~isempty(idx)
            CTD_data.raw_data_plusox.Oxygen(i) = dummy.OXYGEN(idx);
            CTD_data.raw_data_plusox.Oxygen_qc(i) = dummy.OXYGEN_QC(idx);
            CTD_data.raw_data_plusox.psal(i) = dummy.SALINITY(idx);
            CTD_data.raw_data_plusox.psal_qc(i) = dummy.SALINITY_QC(idx);
            CTD_data.raw_data_plusox.temp(i) = dummy.TEMPERATURE(idx);
            CTD_data.raw_data_plusox.temp_qc(i) = dummy.TEMPERATURE_QC(idx);
        else
            CTD_data.raw_data_plusox.Oxygen(i) = NaN;
            CTD_data.raw_data_plusox.Oxygen_qc(i) = NaN;
            CTD_data.raw_data_plusox.psal(i) = NaN;
            CTD_data.raw_data_plusox.psal_qc(i) = NaN;
            CTD_data.raw_data_plusox.temp(i) = NaN;
            CTD_data.raw_data_plusox.temp_qc(i) = NaN;
        end
    end
end

CTD_data.raw_data_plusox.Oxygen(CTD_data.raw_data_plusox.Oxygen==0)=NaN;
% CTD_data.raw_data_plusox.Oxygen_qc(CTD_data.raw_data_plusox.Oxygen_qc==0)=NaN;
CTD_data.raw_data_plusox.psal(CTD_data.raw_data_plusox.psal==0)=NaN;
CTD_data.raw_data_plusox.temp(CTD_data.raw_data_plusox.temp==0)=NaN;
% CTD_data.raw_data_plusox.psal_qc(CTD_data.raw_data_plusox.psal_qc==0)=NaN;
% CTD_data.raw_data_plusox.temp_qc(CTD_data.raw_data_plusox.temp_qc==0)=NaN;

CTD_data.raw_data_plusox(1,:)=[];

Coor_ox_bottleCTD = [CTD_data.raw_data_plusox.Longitude, CTD_data.raw_data_plusox.Latitude, CTD_data.raw_data_plusox.Depth];
Meas_ox_bottleCTD = [CTD_data.raw_data_plusox.psal, CTD_data.raw_data_plusox.Oxygen, CTD_data.raw_data_plusox.temp];
% define the measurement IDs (salinity, oxygen, temperature)
MeasID = [1 6 7];


CTD_data.raw_data_plusox.pH_LIPHR = LIPHR(Coor_ox_bottleCTD, Meas_ox_bottleCTD, MeasID,'OAAdjustTF',false,'PerKgSwTF',false);%oxygen sensor data in umol/L
CTD_data.raw_data_plusox.pH_LIPHR_DIC_Alk = LIPHR(Coor_ox_bottleCTD, Meas_ox_bottleCTD, MeasID,'OAAdjustTF',false,'PerKgSwTF',false,'pHCalcTF',true);%oxygen sensor data in umol/L, and turned on pHCalcTF to get a better alignment with pH calculated from DIC and alk

estdate = years(CTD_data.raw_data_plusox.date - datetime(2000,1,1)) + 2000;
CTD_data.raw_data_plusox.pH_LIPHR_OA = LIPHR(Coor_ox_bottleCTD, Meas_ox_bottleCTD, MeasID, 'EstDates', estdate,'OAAdjustTF', true,'PerKgSwTF',false);%oxygen sensor data in umol/L
CTD_data.raw_data_plusox.pH_LIPHR_OA_DIC_Alk = LIPHR(Coor_ox_bottleCTD, Meas_ox_bottleCTD, MeasID, 'EstDates', estdate,'OAAdjustTF', true,'PerKgSwTF',false);%oxygen sensor data in umol/L, and turned on pHCalcTF to get a better alignment with pH calculated from DIC and alk

% the CTD profiles for 4th and 8th March 2018 seem erroneous, especially
% the oxygen trace, I will take them out

CTD_data.raw_data_plusox.pH_LIPHR(datenum(CTD_data.raw_data_plusox.date)>datenum('04-Mar-2018','dd-mmm-yyyy')&datenum(CTD_data.raw_data_plusox.date)<datenum('09-Mar-2018','dd-mmm-yyyy'),:)=NaN;
CTD_data.raw_data_plusox.pH_LIPHR_DIC_Alk(datenum(CTD_data.raw_data_plusox.date)>datenum('04-Mar-2018','dd-mmm-yyyy')&datenum(CTD_data.raw_data_plusox.date)<datenum('09-Mar-2018','dd-mmm-yyyy'),:)=NaN;
CTD_data.raw_data_plusox.pH_LIPHR_OA(datenum(CTD_data.raw_data_plusox.date)>datenum('04-Mar-2018','dd-mmm-yyyy')&datenum(CTD_data.raw_data_plusox.date)<datenum('09-Mar-2018','dd-mmm-yyyy'),:)=NaN;
CTD_data.raw_data_plusox.pH_LIPHR_OA_DIC_Alk(datenum(CTD_data.raw_data_plusox.date)>datenum('04-Mar-2018','dd-mmm-yyyy')&datenum(CTD_data.raw_data_plusox.date)<datenum('09-Mar-2018','dd-mmm-yyyy'),:)=NaN;


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% all float pH vs CTD casts at SOTS
% load('CTD_data.mat')


figure()
for i = 1:6
    subplot(3,2,i)
    hold on
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==i), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==i),'.b','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR(month(CTD_data_ox.START_TIME)==i), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i), '.k','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR_DIC_Alk(month(CTD_data_ox.START_TIME)==i), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i), '^k','MarkerSize',2)
    plot(CTD_data_ox.pH_LIPHR_OA(month(CTD_data_ox.START_TIME)==i), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i), '.c','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR_OA_DIC_Alk(month(CTD_data_ox.START_TIME)==i), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i), '^c','MarkerSize',2)
   
    if i ==1
        xlabel('pH total scale - Jan')
    elseif i ==2
        xlabel('pH total scale - Feb')
    elseif i==3
        xlabel('pH total scale - Mar')
    elseif i==4
        xlabel('pH total scale - Apr')
    elseif i==5
        xlabel('pH total scale - May')
    elseif i==6
        xlabel('pH total scale - June')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
    %legend('blue CTD','black CTD LIPHR','cyan CTD LIPHR OA','Location','best')

end

figure()
for i = 1:6
    subplot(3,2,i)
    hold on
    plot(CTD_data.raw_data.pH(month(CTD_data.raw_data.date)==(i+6)), CTD_data.raw_data.Depth(month(CTD_data.raw_data.date)==(i+6)),'.b','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR(month(CTD_data_ox.START_TIME)==i+6), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i+6), '.k','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR_DIC_Alk(month(CTD_data_ox.START_TIME)==i+6), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i+6), '^k','MarkerSize',3)
    plot(CTD_data_ox.pH_LIPHR_OA(month(CTD_data_ox.START_TIME)==i+6), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i+6), '.c','MarkerSize',10)
    plot(CTD_data_ox.pH_LIPHR_OA_DIC_Alk(month(CTD_data_ox.START_TIME)==i+6), CTD_data_ox.PRESSURE(month(CTD_data_ox.START_TIME)==i+6), '^c','MarkerSize',3)
 
    if i ==1
        xlabel('pH total scale - July')
    elseif i ==2
        xlabel('pH total scale - Aug')
    elseif i==3
        xlabel('pH total scale - Sep')
    elseif i==4
        xlabel('pH total scale - Oct')
    elseif i==5
        xlabel('pH total scale - Nov')
    elseif i==6
        xlabel('pH total scale - Dec')
    end
    ylabel('Depth dbar')
%     xlim([7.9 8.2])
    ylim([0 2000])
    set(gca, 'YDir','reverse')
%     legend('blue CTD','black CTD LIPHR','cyan CTD LIPHR OA','Location','best')
end


%%%%%%%%% now only LIPHR calculated pH at depths where DIC and Alk were
%%%%%%%%% measured

n=1;
figure()
for i = 1:12
    p1 = plot(CTD_data.raw_data_plusox.pH(month(CTD_data.raw_data_plusox.date)==i & ~isnan(CTD_data.raw_data_plusox.pH_LIPHR)), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i & ~isnan(CTD_data.raw_data_plusox.pH_LIPHR)),'.b','MarkerSize',10); % pH calculated from bottle DIC and Alk
    
    if ~isempty(p1)
        subplot(3,2,n)
        hold on
        plot(CTD_data.raw_data_plusox.pH(month(CTD_data.raw_data_plusox.date)==i & ~isnan(CTD_data.raw_data_plusox.pH_LIPHR)), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i & ~isnan(CTD_data.raw_data_plusox.pH_LIPHR)),'.b','MarkerSize',10) % pH calculated from bottle DIC and Alk
        plot(CTD_data.raw_data_plusox.pH_LIPHR_DIC_Alk(month(CTD_data.raw_data_plusox.date)==i), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i),'.k','MarkerSize',10) % pH calculated from LIPHR based on CTD sensor oxygen, temp and psal, turned on function to align with DIC/Alk calc
        plot(CTD_data.raw_data_plusox.pH_LIPHR_OA_DIC_Alk(month(CTD_data.raw_data_plusox.date)==i), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i),'.c','MarkerSize',10) % pH calculated from LIPHR based on CTD sensor oxygen, temp and psal, with OA correction, turned on function to align with DIC/Alk calc
        xlabel(['pH total scale - Month ', num2str(i)])
        ylabel('Depth dbar')
        set(gca, 'YDir','reverse')
        legend('\color{blue}CTD','\color{black}CTD LIPHR','\color{cyan}CTD LIPHR OA','Location','best')
        ylim([0 2000])
        n=n+1;
    else
        n=n;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% now the difference between DIC / Alk calculated pH and LIPHR estimated
%%% pH based on sensor O2, psal and temp

delta_LIPHR = CTD_data.raw_data_plusox.pH - CTD_data.raw_data_plusox.pH_LIPHR_DIC_Alk; 
% delta_LIPHR(isnan(CTD_data.raw_data_plusox.pH) & isnan(CTD_data.raw_data_plusox.pH_LIPHR_DIC_Alk))=[];
% delta_LIPHR(415:438)=NaN;

delta_LIPHR_OA = CTD_data.raw_data_plusox.pH - CTD_data.raw_data_plusox.pH_LIPHR_OA_DIC_Alk; 
% delta_LIPHR_OA(isnan(CTD_data.raw_data_plusox.pH) & isnan(CTD_data.raw_data_plusox.pH_LIPHR_OA_DIC_Alk))=[];
% delta_LIPHR_OA(415:438)=NaN;


n=1;
figure()
for i = 1:12
    p1 = plot(delta_LIPHR(month(CTD_data.raw_data_plusox.date)==i), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i),'.b','MarkerSize',10); % pH calculated from bottle DIC and Alk - LIPHR estimate
    
    if ~isempty(p1)
        subplot(3,2,n)
        hold on
        plot(delta_LIPHR(month(CTD_data.raw_data_plusox.date)==i), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i),'.b','MarkerSize',10); % pH calculated from bottle DIC and Alk - LIPHR estimate
        plot(delta_LIPHR_OA(month(CTD_data.raw_data_plusox.date)==i), CTD_data.raw_data_plusox.Depth(month(CTD_data.raw_data_plusox.date)==i),'.k','MarkerSize',10); % pH calculated from bottle DIC and Alk - LIPHR OA estimate
        xlabel(['delta pH - Month ', num2str(i)])
        ylabel('Depth dbar')
        set(gca, 'YDir','reverse')
        legend('\color{blue}bottle - LIPHR','\color{black}bottle - LIPHR OA','Location','best')
        ylim([0 2000])
%         xlim([-0.08 0.03])
        n=n+1;
    else
        n=n;
    end
    
end

LIPHR_error = max(abs(delta_LIPHR(delta_LIPHR<=0.03))) % to take out those suspicious March casts
LIPHR_OA_error = max(abs(delta_LIPHR_OA(delta_LIPHR_OA<=0.03))) % to take out those suspicious March casts



%%% find the pH at the deep reference depth (1500m)
ref_depth_D = 1500;
%%% set the tolerance to the reference depth
ref_depth_t_D = 100;
%%% find the pH at the deep reference depth (950m)
ref_depth_S = 950;
%%% set the tolerance to the reference depth
ref_depth_t_S = 100;
%%% find the index of each cast profile where the depth is closest to the
%%% reference depth
%%% first, subset the whole table to only include one cast at a time
n=1;
for i = 1:height(CTD_data.raw_data_plusox)
    l = strcmp(CTD_data.raw_data_plusox.Sample(i),CTD_data.raw_data_plusox.Sample(:)) & CTD_data.raw_data_plusox.CTD(i)==CTD_data.raw_data_plusox.CTD(:);
    dummy = CTD_data.raw_data_plusox(l,:);
    if sum(dummy.Depth>=(ref_depth_D-ref_depth_t_D))>0
        %%% find the pH at the reference depths, 1500m
        pres_ref_id_D = dummy.Depth(:)>=(ref_depth_D-ref_depth_t_D) & dummy.Depth(:)<=(ref_depth_D+ref_depth_t_D); 
        CTD_data.pH_1500(n) = mean(dummy.pH(pres_ref_id_D),'omitnan'); % pH at 1500m, calculated from DIC and Alk
        CTD_data.pH_1500_LIPHR(n) = mean(dummy.pH_LIPHR_DIC_Alk(pres_ref_id_D),'omitnan'); % LIPHR estimate at 1500m, with DIC/Alk adjustment
        CTD_data.pH_1500_LIPHR_OA(n) = mean(dummy.pH_LIPHR_OA_DIC_Alk(pres_ref_id_D),'omitnan'); % LIPHR estimate at 1500m, with OA adjustment, with DIC/Alk adjustment
        %%% find the pH at the reference depths, 950m
        pres_ref_id_S = dummy.Depth(:)>=(ref_depth_S-ref_depth_t_S) & dummy.Depth(:)<=(ref_depth_S+ref_depth_t_S); 
        CTD_data.pH_950(n) = mean(dummy.pH(pres_ref_id_S),'omitnan'); % pH at 950m, calculated from DIC and Alk
        CTD_data.pH_950_LIPHR(n) = mean(dummy.pH_LIPHR_DIC_Alk(pres_ref_id_S),'omitnan'); % LIPHR estimate at 950m, with DIC/Alk adjustment
        CTD_data.pH_950_LIPHR_OA(n) = mean(dummy.pH_LIPHR_OA_DIC_Alk(pres_ref_id_S),'omitnan'); % LIPHR estimate at 950m, with OA adjustment, with DIC/Alk adjustment
        n=n+1;
    else
        n=n;
    end
end


figure()
plot(CTD_data.pH_1500,CTD_data.pH_950,'o')
xlabel('pH at 1500m')
ylabel('pH at 950m')
hold on
plot(CTD_data.pH_1500_LIPHR,CTD_data.pH_950_LIPHR,'.r')
plot(CTD_data.pH_1500_LIPHR_OA,CTD_data.pH_950_LIPHR_OA,'.k')
legend('\color{blue} DIC/Alk calculated pH','\color{red}LIPHR calculated pH','\color{black}LIPHR OA calculated pH','Location','best')


%%%%%%%%%%%%%%%%%%%%%%
%%% compare mooring based estimated pH to float pH
%%% I will use mooring psal and Elizabeth's linear regression plus mooring
%%% pCO2 to calculate pH with COSYS

load('mooring_data.mat')
load('SOTS_float_data_TEMP.mat')
load('WOA_nut_SOTS.mat')
nut_data = WOA_nut_SOTS;

 % calculate Alkalinity based on the linear regression for the top 100m
% from the RAS alk and dic QC report
% alk (umol/kg) = 39.23* salinity + 937.3
Alk_ES_mooring  = 39.23*mooring_data.xCO2_PSAL+ 937.3;
for i = 1:length(mooring_data.pCO2_sw)
    [DATA_L_D, HEADERS, NICEHEADERS] = CO2SYS(mooring_data.pCO2_sw(i),Alk_ES_mooring(i),4,1,mooring_data.xCO2_PSAL(i),mooring_data.xCO2_SST(i),mooring_data.xCO2_SST(i),mooring_data.xCO2_pres(i),mooring_data.xCO2_pres(i),nut_data.silicate.annual_mean,nut_data.phosphate.annual_mean,2,0,1,10,1,2,2);
    mooring_data.mooring_pH(i) = DATA_L_D(43);
end
mooring_data.mooring_pH(mooring_data.mooring_pH==-999)=NaN;

% I also want to see what the LIPHR algorithm would predict from mooring
% potential temperature and salinity (in the absence of O2)

mooring_data.theta = sw_ptmp(mooring_data.xCO2_PSAL,mooring_data.xCO2_SST,mooring_data.xCO2_pres,10);
Coor_mooring = [mooring_data.xCO2_lon, mooring_data.xCO2_lat, mooring_data.xCO2_pres];
Meas_mooring = [mooring_data.xCO2_PSAL, mooring_data.theta];
Meas_mooring_2 = [mooring_data.xCO2_PSAL, mooring_data.xCO2_SST];
% define the measurement IDs (salinity, oxygen, temperature)
MeasID = [1 2];
MeasID_2 = [1 7];

mooring_data.pH_LIPHR = LIPHR(Coor_mooring, Meas_mooring, MeasID,'OAAdjustTF',false);

estdate = years(mooring_data.xCO2_time - datetime(2000,1,1)) + 2000;
mooring_data.pH_LIPHR_OA = LIPHR(Coor_mooring, Meas_mooring, MeasID, 'EstDates', estdate,'OAAdjustTF', true); % including OA estimate
mooring_data.pH_LIPHR_OA_2 = LIPHR(Coor_mooring, Meas_mooring_2, MeasID_2, 'EstDates', estdate,'OAAdjustTF', true); % including OA estimate


%%%%% now let's plot mooring pH calculated from pCO2 and Alk (based on ES salinity) with CO2SYS 

figure()
plot(mooring_data.xCO2_time,mooring_data.mooring_pH,'ok') % this is pH calculated from pCO2 and Alk (based on ES salinity) with CO2SYS 
hold on
for i = 1:length(SOTS_float_data_TEMP.time)
    plot(SOTS_float_data_TEMP.time(i),SOTS_float_data_TEMP.pH_LD_corr(SOTS_float_data_TEMP.pres(:,i)<20,i),'ob','MarkerSize',2)
end
for i = 1:length(SOTS_float_data_TEMP.time)
    plot(SOTS_float_data_TEMP.time(i),SOTS_float_data_TEMP.pH_LS_corr(SOTS_float_data_TEMP.pres(:,i)<20,i),'.r')
end
plot(SOTS_float_data_TEMP.time, SOTS_float_data_TEMP.pH_L_S_20_corr,'.g')
for i = 1:length(time_float)
    plot(time_float(i),pH_adjusted(pres(:,i)<20,i),'*c','MarkerSize',2)
end
plot(mooring_data.xCO2_time,mooring_data.pH_LIPHR_OA,'om') % this is pH calculated from salinity and potential temperature with LIPHR
plot(mooring_data.xCO2_time,mooring_data.pH_LIPHR_OA,'og') % this is pH calculated from salinity and potential temperature with LIPHR

xlim([min(SOTS_float_data_TEMP.time), max(SOTS_float_data_TEMP.time)])
legend('\color{black} mooring pH from pCO2 and Alk','\color{blue} float LD','\color{red} float LS','\color{green} float LS top 20m ave','\color{cyan} float pH QCed', '\color{green} mooring pH from LIPHR OA') 



figure()
plot(mooring_data.xCO2_time,mooring_data.mooring_pH,'ok')
hold on
plot(SOTS_float_data_TEMP.time, SOTS_float_data_TEMP.pH_L_S_20_corr,'.-g')
plot(SOTS_float_data_TEMP.time, SOTS_float_data_TEMP.pH_L_D_20_corr,'.-b')
xlim([min(SOTS_float_data_TEMP.time), max(SOTS_float_data_TEMP.time)])
legend('\color{black} mooring pH','\color{green} float LS top 20m ave','\color{blue} float LD top 20m ave') 


%%%%%%%%%%%
%%% assuming that float and mooring sample the same watermasses in time,
%%% does the difference between LS/LD pH and mooring pH correlate with any
%%% parameter, such as difference in lat/lon?

for i = 1:length(SOTS_float_data_TEMP.pH_L_S_20_corr)
    idx = knnsearch(datenum(mooring_data.xCO2_time(:)),datenum(SOTS_float_data_TEMP.time(i)));
    mooring_data.mooring_pH_float_overlap(i) = mooring_data.mooring_pH(idx);
    mooring_data.mooring_pH_float_overlap_time(i) = mooring_data.xCO2_time(idx);
    mooring_data.mooring_PSAL_float_overlap(i) = mooring_data.xCO2_PSAL(idx);
    mooring_data.mooring_TEMP_float_overlap(i) = mooring_data.xCO2_SST(idx);
end

figure()
% plot(mooring_data.xCO2_time,mooring_data.mooring_pH,'.k')
hold on
plot(mooring_data.mooring_pH_float_overlap_time,mooring_data.mooring_pH_float_overlap,'ob')
plot(SOTS_float_data_TEMP.time, SOTS_float_data_TEMP.pH_L_S_20_corr,'og')
plot(SOTS_float_data_TEMP.time, SOTS_float_data_TEMP.pH_L_D_20_corr,'om')
xlim([min(SOTS_float_data_TEMP.time), max(SOTS_float_data_TEMP.time)])
legend('\color{blue} mooring pH','\color{green} float LS top 20m ave','\color{magenta} float LD top 20m ave') 

mooring_data.pH_diff_LS = mooring_data.mooring_pH_float_overlap - SOTS_float_data_TEMP.pH_L_S_20_corr;
mooring_data.pH_diff_LD = mooring_data.mooring_pH_float_overlap - SOTS_float_data_TEMP.pH_L_D_20_corr;

figure()
hold on
plot(SOTS_float_data_TEMP.time, mooring_data.pH_diff_LS, 'ob')
plot(SOTS_float_data_TEMP.time, mooring_data.pH_diff_LD, 'om')
% plot(SOTS_float_data_TEMP.time, mooring_data.pH_diff_LS-mooring_data.pH_diff_LD,'og')
yline([0])
legend('\color{blue} pH diff LS','\color{magenta} pH diff LD','\color{green} pH diff between the two') 

figure()
plot(mooring_data.pH_diff_LS,SOTS_float_data_TEMP.mooring_lat-SOTS_float_data_TEMP.lat,'.k')


mdl_LS_lat = fitlm(mooring_data.pH_diff_LS,SOTS_float_data_TEMP.mooring_lat-SOTS_float_data_TEMP.lat)
mdl_LD_lat = fitlm(mooring_data.pH_diff_LD,SOTS_float_data_TEMP.mooring_lat-SOTS_float_data_TEMP.lat)

figure()
plot(mdl_LS_lat)
xlabel('pH mooring - float pH_LS')
ylabel('mooring lat - float lat')

figure()
plot(mdl_LD_lat)
xlabel('pH mooring - float pH_LD')
ylabel('mooring lat - float lat')



mdl_LS_lon = fitlm(mooring_data.pH_diff_LS,SOTS_float_data_TEMP.mooring_lon-SOTS_float_data_TEMP.lon)
mdl_LD_lon = fitlm(mooring_data.pH_diff_LD,SOTS_float_data_TEMP.mooring_lon-SOTS_float_data_TEMP.lon)

figure()
plot(mdl_LS_lon)
xlabel('pH mooring - float pH_LS')
ylabel('mooring lon - float lon')

figure()
plot(mdl_LD_lon)
xlabel('pH mooring - float pH_LD')
ylabel('mooring lon - float lon')



%%%% compare T & S of mooring with that of float, are they sampling the
%%%% same water masses?

% all mooring data points
figure()
plot(mooring_data.xCO2_PSAL, mooring_data.xCO2_SST,'ob')
hold on
plot(SOTS_float_data_TEMP.PSAL_20, SOTS_float_data_TEMP.TEMP_20, '.r')
xlabel('PSAL')
ylabel('TEMP')

% only the mooring data points that coincide in TIME with the float
% profiles
figure()
plot(mooring_data.mooring_PSAL_float_overlap, mooring_data.mooring_TEMP_float_overlap,'ob')
hold on
plot(SOTS_float_data_TEMP.PSAL_20, SOTS_float_data_TEMP.TEMP_20, '.r')
xlabel('PSAL')
ylabel('TEMP')

% only the mooring data points that coincide in TIME with the float
% profiles and a colorbar to indicate what time those points are
c=datenum(mooring_data.mooring_pH_float_overlap_time);
cnum=[datenum('01-Jan-2021','dd-mmm-yyyy'), datenum('01-Mar-2021','dd-mmm-yyyy'),... 
    datenum('01-Jun-2021','dd-mmm-yyyy'), datenum('01-Sep-2021','dd-mmm-yyyy'),...
    datenum('01-Dec-2021','dd-mmm-yyyy')];
clabel={'Jan', 'Mar', 'Jun', 'Sep', 'Dec'};
d=datenum(SOTS_float_data_TEMP.time);

figure()
yyaxis left
scatter(mooring_data.mooring_PSAL_float_overlap, mooring_data.mooring_TEMP_float_overlap,20,c,'filled')
colorbar('Ticks',cnum,'TickLabels',clabel,'Location','west')
ylim([8 14])

hold on
yyaxis right
scatter(SOTS_float_data_TEMP.PSAL_20, SOTS_float_data_TEMP.TEMP_20,30,d,'p','filled')
colorbar('Ticks',cnum,'TickLabels',clabel,'Location','east')
ylim([8 14])
xlabel('PSAL')
ylabel('TEMP')



%%%% does diff pH correlate with diff T&S?
% what if I limit the data points to compare to those where mooring data
% and float data are within 1 or 2 days and with temp and psal limits that
% equate to those set for well mixed waters (mixed layer depth
% definitions), to make sure we only compare data from the same water
% masses. Temp threshold = 0.3C (temp QC report), psal threshold = 0.03psu
% (salinity QC report)

% time threshold in days
d = 2;
% temp threshold in C
t = 0.3;
% psal threshold in psu
p = 0.03;

for i = 1:size(SOTS_float_data_TEMP.pH,2)
    %subset by time threshold
    temporal_sub = abs(datenum(mooring_data.xCO2_time)-datenum(SOTS_float_data_TEMP.time(i)))<=d;
    %subset by temp threshold
    temp_sub = abs(mooring_data.xCO2_SST-SOTS_float_data_TEMP.TEMP_20(i))<=t;
    %subset by psal threshold
    psal_sub = abs(mooring_data.xCO2_PSAL-SOTS_float_data_TEMP.PSAL_20(i))<=p;
    %now take the mean of the mooring pH values where all there indices = 1
    mooring_pH_subset(i) = mean(mooring_data.mooring_pH(temporal_sub==1 & temp_sub ==1 & psal_sub==1),'omitnan');
    mooring_flux_subset(i) = mean(mooring_data.flux_pCO2(temporal_sub==1 & temp_sub ==1 & psal_sub==1),'omitnan');
    
end

figure()
plot(SOTS_float_data_TEMP.pH_L_S_20_corr,mooring_pH_subset,'*k')
xlabel('float pH')
ylabel('mooring pH')
hold on
plot(SOTS_float_data_TEMP.pH_L_D_20_corr,mooring_pH_subset,'+r')
ylim([8.03 8.11])
plot([8.03 8.11],[8.03 8.11],'-k')
legend('LIPHR Shallow', 'LIPHR Deep')

res_LS_pH = sqrt(sum((SOTS_float_data_TEMP.pH_L_S_20_corr-mooring_pH_subset).^2,'omitnan'));
res_LD_pH = sqrt(sum((SOTS_float_data_TEMP.pH_L_D_20_corr-mooring_pH_subset).^2,'omitnan'));


figure()
plot(SOTS_float_data_TEMP.flux_L_S_corr,mooring_flux_subset,'*k')
xlabel('float air sea flux')
ylabel('mooring air sea flux')
hold on
plot(SOTS_float_data_TEMP.flux_L_D_corr,mooring_flux_subset,'+r')
ylim([-30 5])
plot([-30 5],[-30 5],'-k')
legend('LIPHR Shallow', 'LIPHR Deep')

res_LS_flux = sqrt(sum((SOTS_float_data_TEMP.flux_L_S_corr-mooring_flux_subset).^2,'omitnan'));
res_LD_flux = sqrt(sum((SOTS_float_data_TEMP.flux_L_D_corr-mooring_flux_subset).^2,'omitnan'));
