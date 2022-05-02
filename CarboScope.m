%%%%%%%%%%%%%%% a script to read the CarboScope data
clear all
load('S55_float_data.mat')

file ='C:\Users\cawynn\cloudstor\Air sea flux manuscript\CarboScope\oc_v2021_daily.nc';

CS.time = seconds(ncread(file, 'mtime'))+ datetime(2000,1,1);
CS.flux = ncread(file,'co2flux_ocean'); % PgC/yr
CS.dxyp = ncread(file,'dxyp'); % divide by this to convert PgC/yr from grid cell to m-2
CS.lat = ncread(file,'lat');
CS.lon = ncread(file,'lon');


for i = 1:size(S55_float_data.lat,2)
%for i = 56   
    % find the CarboScope file index that is closest in time and space to the float
    idx_T = knnsearch(datenum(CS.time(:)),datenum(S55_float_data.time(i)));
    % delta_time = (EraTime(idx_T)-time_float(i))/24;
    idx_lat = knnsearch(CS.lat(:),S55_float_data.lat(i));
    idx_lon = knnsearch(CS.lon(:),S55_float_data.lon(i));
    
    CS.S55_dxyp(i) = CS.dxyp(idx_lon,idx_lat);
    CS.S55_flux_PgC(i) = CS.flux(idx_lon,idx_lat,idx_T)/CS.S55_dxyp(i); % PgC/yrm2
    CS.S55_flux_molC(i) = CS.S55_flux_PgC(i)/12*10^15; % mol C / y / m2 
    CS.S55_time(i) = CS.time(idx_T);
end



CS.clim_mo = month(CS.time(15706:23376));
CS.clim=CS.flux(128:138,16:19,15706:23376)./CS.dxyp(128:138,16:19);%PgC/y/m2

% sort the climatology by month and then take the monthly means
[a_sorted a_order]= sort(CS.clim_mo);
CS.clim_sorted = CS.clim(:,:,a_order);

for i=1:12
    CS.clim_flux(i) = mean(CS.clim_sorted(:,:,a_sorted==i),...
        'all','omitnan')/12*10^15; % mol C y-1 m-2
    
    CS.clim_flux_SD(i) = std(CS.clim_sorted(:,:,a_sorted==i),0,...
        [1 2 3 4],'omitnan')/12*10^15; % mol C y-1 m-2
end

clearvars -except CS
save('CS.mat')

load('S55_float_data.mat')

figure()
plot(S55_float_data.time, S55_float_data.flux_L_S_corr,'*b')
hold on
plot(CS.S55_time, CS.S55_flux_molC,'or')
xlabel('time')
ylabel('air sea flux mol C m^-^2 y^-^1')
legend('S55 float 2021', 'CarboScope')

figure()
plot([1:12],CS.clim_flux,'.k')
xlabel('time')
ylabel('air sea flux mol C m^-^2 y^-^1')
legend('CarboScope climatology 2000-2020')
       

%%% plot the different air sea fluxes vs S55 float as monthly means
load('S55_comp_float_data.mat')
load('S55_float_data.mat')
addpath 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts'

figure()
plot([1:12],(S55_comp_float_data.y2016.flux_pH_adj_mo_ave/1000)*365,'*r')
hold on
plot([1:12],(S55_comp_float_data.y2016.flux_LD_mo_ave/1000)*365,'*b')
plot([1:12],(S55_comp_float_data.y2016.flux_LS_mo_ave/1000)*365,'*g')
plot([1:12],(S55_comp_float_data.y2017.flux_pH_adj_mo_ave/1000)*365,'or')
plot([1:12],(S55_comp_float_data.y2017.flux_LD_mo_ave/1000)*365,'ob')
plot([1:12],(S55_comp_float_data.y2017.flux_LS_mo_ave/1000)*365,'og')
plot(S55_comp_float_data.y2018.mo_ave_month,(S55_comp_float_data.y2018.flux_pH_adj_mo_ave/1000)*365,'.r')
plot(S55_comp_float_data.y2018.mo_ave_month,(S55_comp_float_data.y2018.flux_LD_mo_ave/1000)*365,'.b')
plot(S55_comp_float_data.y2018.mo_ave_month,(S55_comp_float_data.y2018.flux_LS_mo_ave/1000)*365,'.g')
plot([1:12],(S55_float_data.flux_LD_mo_ave/1000)*365,'oc','MarkerFaceColor','c','MarkerSize',6)
plot([1:12],(S55_float_data.flux_LS_mo_ave/1000)*365,'om','MarkerFaceColor','m','MarkerSize',6)
% plot([1:12],CS.clim_flux,'.k','MarkerSize',10)
shadedErrorBar([1:12],CS.clim_flux,...
        CS.clim_flux_SD,'lineprops','-.k','patchSaturation',[0.03])
hold off
ylabel('air sea flux mol m^-^2 y^-^1')
legend('red * pH adj 2016', 'blue * pH LD 2016', 'green * pH LS 2016',...
'red o pH adj 2017', 'blue o pH LD 2017', 'green o pH LS 2017',...
'red . pH adj 2018', 'blue . pH LD 2018', 'green . pH LS 2018',...
'cyan + pH LD S55', 'magenta + pH LS S55','CS climatology 2000-2020',...
'Orientation','vertical','Location','westoutside')
xticks([0:13])
xticklabels({'', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug',...
    'Sep', 'Oct', 'Nov', 'Dec',''}) 
xlim([0 13])
yline(0)