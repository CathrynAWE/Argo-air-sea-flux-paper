
file='woa18_all_n00_01.nc.save';
var_to_plot = 'n_dd';

figure(1); clf

    
%file='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/nitrate/all/1.00/woa18_all_n00_01.nc';
%file='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/silicate/all/1.00/woa18_all_i00_01.nc';

% get the file info
f_info = ncinfo(file);

% get the lat and lon grid points
lat_grid=ncread(file, 'lat');
lon_grid=ncread(file,'lon');

lat=[-60 -30];
lon=[120 160];

% find the index of the lat/lon grid points
lat_idx(1) = find(lat_grid>lat(1),1);
lon_idx(1) = find(lon_grid>lon(1),1);
lat_idx(2) = find(lat_grid>lat(2),1);
lon_idx(2) = find(lon_grid>lon(2),1);

lat_sots=-47;
lon_sots=142;

lat_idx_sots = find(lat_grid>lat_sots,1);
lon_idx_sots = find(lon_grid>lon_sots,1);

depth=0;

depth_grid = ncread(file, 'depth');
depth_idx = find(depth_grid>=depth,1);

t = ncread(file, 'time');
clear data;

var_name = ncreadatt(file, var_to_plot, 'standard_name');
var_units = ncreadatt(file, var_to_plot, 'units');

for i=1:numel(f_info.Variables)
    if isfield(f_info.Variables(i), 'Dimensions')
        if numel(f_info.Variables(i).Dimensions) >= 4
            name = f_info.Variables(i).Name;

            std_name_idx = find(strcmp({f_info.Variables(i).Attributes.Name}, 'standard_name')==1);
            if numel(std_name_idx) == 0
                std_name_idx = find(strcmp({f_info.Variables(i).Attributes.Name}, 'long_name')==1);
            end

            %disp(f_info.Variables(i).Attributes(std_name_idx).Value)

            % read and save data

            %data.(name) = ncread(file, name, [lon_idx lat_idx depth_idx 1], [2 2 1 1]);            
            data.(name) = ncread(file, name);
        end
    end
end

xy = permute(data.(var_to_plot)((lon_idx(1):lon_idx(2)),(lat_idx(1):lat_idx(2)),1),[2 1]);

imagesc(lon_grid(lon_idx(1):lon_idx(2)), lat_grid(lat_idx(1):lat_idx(2)), xy(:,:,1), 'AlphaData',(~isnan(xy(:,:,1))))
title(var_name, 'Interpreter', 'None', 'FontWeight', 'normal');
colorbar
set(gca,'YDir','normal')
caxis([0 15]);

