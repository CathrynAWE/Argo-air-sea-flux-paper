%%%% for the CO2SYS functions I need nutrients, here I will extract them
%%%% for the float profiles from WOA 2018
clear all

% float_ID = '5906623' % SOTS float
float_ID = '5906624' %55S float

if float_ID == '5906623'
    load('SOTS_float_data.mat')
    data=SOTS_float_data;
elseif float_ID == '5906624'
    load('S55_float_data.mat')
    data=S55_float_data;
end


lat_range_float = [min(data.lat) max(data.lat)];
lon_range_float = [min(data.lon) max(data.lon)];
depth_range = [0 20];

% get the monthly means
variables ={'silicate','phosphate'};
for l=1:length(variables)
    file_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\WOA\' variables{l}];
    cd(file_path)
    files = dir('*.nc');
    
    if strcmp(variables{l},'silicate')==1
        v='i';
    elseif strcmp(variables{l},'phosphate')==1
        v='p';
    end
    
    for i=1:length(files)
        % find the file of the month we are up to
        if i<10
            pat = [v '0' num2str(i)];
        else
            pat = [v num2str(i)];
        end
        
%         if strcmp(variables{l},'silicate')==1
%             if i<10
%                 pat = [v '0' num2str(i)];
%             else
%                 pat = [v num2str(i)];
%             end
%         elseif strcmp(variables{l},'phosphate')==1
%             if i<10
%                 pat = ['p0' num2str(i)];
%             else
%                 pat = ['p' num2str(i)];
%             end
%         end

        for s=1:length(files)
            p(s)=isempty(strfind(files(s).name,pat));
        end

        f_idx = find(p==0);
        
        
        file = [files(f_idx).folder '\' files(f_idx).name];
        
        
        var_name = [v '_an'];
        lat=ncread(file,'lat');
        lon=ncread(file,'lon');
        depth = ncread(file,'depth');
        
        idx_lon_1 = find(lon>lon_range_float(1),1);
        idx_lon_2 = find(lon>lon_range_float(2),1);
        
        idx_lat_1 = find(lat>lat_range_float(1),1);
        idx_lat_2 = find(lat>lat_range_float(2),1);
        
        depth_1 = find(depth>=depth_range(1),1);
        depth_2 = find(depth>depth_range(2),1);
        
        var=ncread(file,var_name); % umol kg-1
        
        WOA_nut.(variables{l}).mo_mean(i) = mean(var([idx_lon_1:idx_lon_2],...
            [idx_lat_1:idx_lat_2],[depth_1:depth_2]),'all','omitnan');

        
    end
end 

% get the annual mean
variables ={'silicate','phosphate'};
for l=1:length(variables)
    file_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\WOA\' variables{l} '\annual_mean'];
    cd(file_path)
    
    if strcmp(variables{l},'silicate')==1
        v='i';
    elseif strcmp(variables{l},'phosphate')==1
        v='p';
    end
    files = dir(['*' v '00_01*.nc']);
    file = [files.folder '\' files.name];
         

    var_name = [v '_an'];
    lat=ncread(file,'lat');
    lon=ncread(file,'lon');
    depth = ncread(file,'depth');

    idx_lon_1 = find(lon>lon_range_float(1),1);
    idx_lon_2 = find(lon>lon_range_float(2),1);

    idx_lat_1 = find(lat>lat_range_float(1),1);
    idx_lat_2 = find(lat>lat_range_float(2),1);

    depth_1 = find(depth>=depth_range(1),1);
    depth_2 = find(depth>depth_range(2),1);

    var=ncread(file,var_name); % umol kg-1

    WOA_nut.(variables{l}).annual_mean = mean(var([idx_lon_1:idx_lon_2],...
        [idx_lat_1:idx_lat_2],[depth_1:depth_2]),'all','omitnan');
        
end

WOA_nut.lat_range = lat_range_float;
WOA_nut.lon_range = lon_range_float;

path='C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper';
cd(path)

if float_ID == '5906623'
    WOA_nut_SOTS = WOA_nut;
    clearvars -except WOA_nut_SOTS
    save('WOA_nut_SOTS.mat');
elseif float_ID == '5906624'
    WOA_nut_S55 = WOA_nut;
    clearvars -except WOA_nut_S55
    save('WOA_nut_S55.mat')
end

        
% 
% figure()
% plot([1:12],WOA_nut_SOTS.silicate.mo_mean,'or')
% hold on
% plot([1:12],WOA_nut_SOTS.phosphate.mo_mean,'ob')
% plot([1:12],WOA_nut_S55.silicate.mo_mean,'*r')
% plot([1:12],WOA_nut_S55.phosphate.mo_mean,'*b')
        
        