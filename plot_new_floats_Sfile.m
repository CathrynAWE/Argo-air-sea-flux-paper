% plot data from new floats
% light colours are early data, dark colours later in time
% CS, 9.10.2019

% now go find the data
% float IDs, choose one: 5906623, 5906624
float_ID = '5906623'
search_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(search_path)

% find the Sprof
%nc_files = files_w_ext(pwd,'nc'); 
nc_files = dir('*.nc')
the_file = [nc_files.folder '\' nc_files.name]


% for k=1:length(nc_files)
%     blah = strfind(nc_files{k},'Sprof');
%     if ~isempty(blah)
%         ind = k;
%     end
% end

% the_file = nc_files{ind};
% the_file = '5905441_Sprof.nc';

% open the file and extract data of interest
dat = netcdf.open(the_file,'NOWRITE'); 

% find how long the profile is
pp = ncread(the_file,'PRES');

lat = ncread(the_file,'LATITUDE');
lon = ncread(the_file,'LONGITUDE');
time = ncread(the_file,'JULD');%+ datetime(1950,1,1);
pres = ncread(the_file,'PRES');
fluo = ncread(the_file,'CHLA');
bbp700 = ncread(the_file,'BBP700');
%bbp470 = ncread(the_file,'BBP470');    % TEMPO floats have it, SOLACE not
CDOM = ncread(the_file,'CDOM');        % SOLACE floats have it, TEMPO not
oxy = ncread(the_file,'DOXY_ADJUSTED');
NO3 = ncread(the_file,'NITRATE');
pH = ncread(the_file,'PH_IN_SITU_TOTAL_ADJUSTED');
pH_error = ncread(the_file,'PH_IN_SITU_TOTAL_ADJUSTED_ERROR');
%pH = ncread(the_file,'PH_IN_SITU_TOTAL');


PAR = ncread(the_file,'DOWNWELLING_PAR');
Ed380 = ncread(the_file,'DOWN_IRRADIANCE380');
Ed412 = ncread(the_file,'DOWN_IRRADIANCE412');
Ed490 = ncread(the_file,'DOWN_IRRADIANCE490');

T = ncread(the_file,'TEMP');
S = ncread(the_file,'PSAL');

netcdf.close(dat);

fs = size(fluo);

% set colour scheme
C = repmat(linspace(0.8,0,fs(1,2)).',1,3);


%% loop through and plot

close all
for i = 1:fs(1,2)
   
    % how many subplots?
    n=6;
    
    % and plot
    fig = figure(1)
    subplot(1,n,1)
    scatter(fluo(:,i), pres(:,i),15,C(i,:),'filled')
    title('Chl')
    set(gca,'YDir','reverse','XAxisLocation','top')
    ylabel('Pressure')
    hold on
    ylim([0 300])
    set(gca,'box','on')
    %xlim([0 5])

    subplot(1,n,2)
    scatter(bbp700(:,i), pres(:,i),15,C(i,:),'filled')
    title('Backscatter at 700 nm')
    set(gca,'YDir','reverse','XAxisLocation','top')
    hold on
    ylim([0 300])
    set(gca,'box','on')

    subplot(1,n,6)
    scatter(pH(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
    title('pH')
    set(gca,'YDir','reverse','XAxisLocation','top')
    %hold on
    set(gca,'box','on')
    ylim([0 2000])
    xlim([7.5 8.5])

    subplot(1,n,4)
    scatter(oxy(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
    title('Oxygen')
    set(gca,'YDir','reverse','XAxisLocation','top')
    hold on
    set(gca,'box','on')
    grid on
    ylim([0 2000])
    
    subplot(1,n,5)
    scatter(NO3(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
    title('Nitrate')
    set(gca,'YDir','reverse','XAxisLocation','top')
    hold on
    set(gca,'box','on')
    ylim([0 2000])
    
    subplot(1,n,3)
    scatter(CDOM(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
    title('CDOM')
    set(gca,'YDir','reverse','XAxisLocation','top')
    hold on
    set(gca,'box','on')
 
%     subplot(1,n,3)
%     %scatter(bbp470(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
%     title('Backscatter at 470 nm')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     hold on
%     set(gca,'box','on')   
%     ylim([0 300])
%     box on

    %pause
    
%     % make an additional figure for the irradiance sensors ============
%     figure(2)
%     subplot(1,4,1)
%     scatter(PAR(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
%     title('PAR')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     hold on
%     set(gca,'box','on')
%     ylim([0 150])
%     
%     subplot(1,4,2)
%     scatter(Ed380(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
%     title('Ed380')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     hold on
%     set(gca,'box','on')
%     ylim([0 150])
%     
%     subplot(1,4,3)
%     scatter(Ed412(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
%     title('Ed412')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     hold on
%     set(gca,'box','on')
%     ylim([0 150])
%     
%     subplot(1,4,4)
%     scatter(Ed490(:,i), pres(:,i),15,C(i,:),'filled') % needs to be scatter because of NaNs
%     title('Ed490')
%     set(gca,'YDir','reverse','XAxisLocation','top')
%     hold on
%     set(gca,'box','on')
%     ylim([0 150])
%     
%     
%     pause
end
    
orient landscape
figure_path = [search_path '\' 'figures'];
cd(figure_path);
if strcmp(float_ID, '5906623') 
    figure_nm = [float_ID '_' 'SOLACE'];
elseif strcmp(float_ID, '5906624')
    figure_nm = [float_ID '_' 'SOLACE'];
else
    figure_nm = float_ID;
end

saveas(fig, figure_nm,'png')

