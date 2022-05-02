% convert float salinity into alkalinity
% now go find the data
% float IDs, choose one: 5906623, 5906624
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'
float_ID = '5906623'
search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% find the Sprof
%nc_files = files_w_ext(pwd,'nc'); 
nc_files = dir('*.nc')
file = [nc_files.folder '\' nc_files.name]


%dat = netcdf.open(the_file,'NOWRITE'); 

% find how long the profile is
pp = ncread(file,'PRES');

lat = ncread(file,'LATITUDE');
lon = ncread(file,'LONGITUDE');
time = ncread(file,'JULD');%+ datetime(1950,1,1);
pres = ncread(file,'PRES');
oxy = ncread(file,'DOXY_ADJUSTED');
NO3 = ncread(file,'NITRATE');
pH = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH_error = ncread(file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
N = ncread(file, 'NITRATE_ADJUSTED');
T = ncread(file,'TEMP');
S = ncread(file,'PSAL');
fs = size(pH);


%oxy_interp = interp1(pp(~isnan(oxy(:,1)),1),oxy(~isnan(oxy(:,1)),1),pp(:,1));

% LIAR equation masks points where there are NANs in any of the
% parameters
figure()
for i = 1:fs(2)
%for i =1    
    msk = ~isnan(oxy(:,i));
    lon_c = zeros(size(oxy(msk,i),1),1);
    lon_c(:) = lon(i,1);
    lat_c = zeros(size(oxy(msk,i),1),1);
    lat_c(:) = lat(i,1);
    Coor = [lon_c lat_c pp(msk,i)];
%     lon_c = zeros(fs(1),1);
%     lon_c(:) = lon(i,1);
%     lat_c = zeros(fs(1),1);
%     lat_c(:) = lat(i,1);
%     Coor = [lon_c lat_c pp(:,i)];
    Meas = [S(msk,i) oxy(msk,i) T(msk,i)];
%     Meas = [S(:,i) oxy(:,i) T(:,i)];
    MeasID = [1 6 7];

    Alk_LIAR = LIAR(Coor, Meas, MeasID);
    
   
    plot(Alk_LIAR, pp(msk,i),'o')
%     plot(Alk_LIAR, pp(:,i),'o')
    set(gca, 'YDir','reverse')
    ylabel('Pressure')
    xlabel('Alkalinity umol/kg - LIAR')
    hold on
end


%netcdf.close(dat);
% after Roden, et al 2016 - TA+N(umol kg-1) = 67+-1 * S + 36+-18 
Alk_R = 67.*S + 36 - N; % NANs in nitrate data need to be taken out
% use mask for plotting
mskAlk = ~isnan(N);

% after Takahashi, et al 2014 - TA + N(uequ kg-1) = 23.76 * S + 1486.1
Alk_T = 23.76 .*S + 1486.1 -N;


% figure()
% plot(Alk_LIAR, pp(msk,1),'o')
% set(gca, 'YDir','reverse')
% ylabel('Pressure')
% xlabel('Alkalinity umol/kg - LIAR')
% 
% figure()
% for i = 1:fs(2)
% %for i =1
%     plot(Alk_R(mskAlk(:,i),i), pp(mskAlk(:,i),i),'o')
%     set(gca, 'YDir','reverse')
%     ylabel('Pressure')
%     xlabel('Alkalinity umol/kg Roden, 2016')
%     hold on
% end
% 
% figure()
% for i = 1:fs(2)
% %for i =1
%     plot(Alk_T(mskAlk(:,i),i), pp(mskAlk(:,i),i),'o')
%     set(gca, 'YDir','reverse')
%     ylabel('Pressure')
%     xlabel('Alkalinity umol/kg Takahashi, 2014')
%     hold on
% end