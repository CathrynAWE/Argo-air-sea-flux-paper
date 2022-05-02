% this script calculates all of SOFS mooring air-sea fluxes

% calculate flux from dfCO2, SST and S in the files


% patch together daily files, with hourly obs for windspeed and flags
% http://thredds.aodn.org.au/thredds/catalog/IMOS/DWM/ASFS/SOFS/Surface_fluxes/Real-time/2020_daily/catalog.html
clear all

addpath ('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper')

path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\complete_SOFS_xCO2';

cd(path)

nc_files = dir('*.nc');

raw_data= [];
raw_data.wsp =[];
raw_data.wsp_qc =[];
raw_data.wsp_time =[];
raw_data.wsp_lat =[];
raw_data.wsp_lon =[];
raw_data.wsp_airtemp =[];
raw_data.wsp_airtemp_qc = [];

%%% put all the seperate windspeed / air temp nc files together
for i = 1:length(nc_files)
    
    file = [nc_files(i).folder '\' nc_files(i).name];

    wsp = ncread(file,'WSPD10M');
    wsp_qc = ncread(file,'WIND_FLAG');
    wsp_time = ncread(file, 'TIME') + datetime(1950,1,1);
    wsp_lat = ncread(file, 'LATITUDE');
    wsp_lon = ncread(file, 'LONGITUDE');
    wsp_airtemp = ncread(file, 'AIRT');
    wsp_airtemp_qc = ncread(file, 'AIRT_FLAG');
    
    raw_data.wsp = [raw_data.wsp; wsp(:)];
    raw_data.wsp_qc = [raw_data.wsp_qc; wsp_qc(:)];
    raw_data.wsp_time = [raw_data.wsp_time; wsp_time(:)];
    raw_data.wsp_lat = [raw_data.wsp_lat; wsp_lat(:)];
    raw_data.wsp_lon = [raw_data.wsp_lon; wsp_lon(:)];
    raw_data.wsp_airtemp = [raw_data.wsp_airtemp; wsp_airtemp(:)];
    raw_data.wsp_airtemp_qc = [raw_data.wsp_airtemp_qc; wsp_airtemp_qc(:)];

end

% download data from
% https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0118546/
% stitch CSV files together, this is QCed data

path ='C:\Users\cawynn\cloudstor\Air sea flux manuscript\complete_SOFS_xCO2\pCO2';
cd(path)
csv_files=dir('*.csv');


%%% put the separate csv files together
pCO2_data=[];

for i=1:length(csv_files)
    
    file= [csv_files(i).folder '\' csv_files(i).name];
    
    data= readtable(file,'ReadVariableNames',false);
    % some files have two end columns for pH data (always empty), others do
    % not
    if size(data,2)>24
        data(:,25:26)=[];
    end
    
    pCO2_data=[pCO2_data; data];
    
end

%%% put the correct variable names in
dummy=readtable([csv_files(1).folder '\' csv_files(1).name]);
names= dummy.Properties.VariableNames;
if size(names,2)>24
    names(25:26)=[];
end
pCO2_data.Properties.VariableNames = names;

%%% date and time of day are in seperate columns, here I put them back
%%% together, also some csv files don't follow the yyyy format, but instead
%%% just list 15, instead of 2015, for example.

date=string(pCO2_data.Date);
date_new = strrep(date,'/00','/20');
hour=string(pCO2_data.Time);
date_str = strcat(date_new, " ", hour);
pCO2_data.DateTime = datetime(date_str, 'InputFormat', 'MM/dd/yyyy HH:mm');

%%%% sort the table by date
mooring_data=[];
mooring_data.CO2= sortrows(pCO2_data,size(pCO2_data,2),'ascend');

%%% turn -999 into NaNs or else we will get strange air sea flux results
mooring_data.CO2.dpCO2(mooring_data.CO2.dpCO2==-999)=NaN;


% now limit the windspeed readings to the number of pCO2 readings

l=[];
for i = 1:length(mooring_data.CO2.DateTime)
    l = knnsearch(datenum(raw_data.wsp_time(:)),datenum(mooring_data.CO2.DateTime(i)));

    mooring_data.CO2.wsp(i) = raw_data.wsp(l);
    mooring_data.CO2.wsp_time(i) = raw_data.wsp_time(l);
    mooring_data.CO2.wsp_lat(i) = raw_data.wsp_lat(l);
    mooring_data.CO2.wsp_lon(i) = raw_data.wsp_lon(l);
    mooring_data.CO2.wsp_airtemp(i) = raw_data.wsp_airtemp(l);
end


%%% there are times when the mooring wsp data was more than 2h (even gaps
%%% of weeks) apart from the mooring pCO2 observation, for those times I 
%%% will fill in with ERA5 data

% read the ERA5 data for the mooring
load('ERA5_SOFS.mat')

% find the gaps
d=abs(datenum(mooring_data.CO2.DateTime)-datenum(mooring_data.CO2.wsp_time))>2/24;

% fill in the gaps from ERA5
fields=fieldnames(ERA5);
n=length(fields);
s=1;
for i=1:length(d)
    if d(i)==1
        for s=1:n
            idx = knnsearch(datenum(getfield(ERA5,fields{s},'Time')),datenum(mooring_data.CO2.DateTime(i)));
            if datenum(getfield(ERA5,fields{s},'Time',{idx}))-datenum(mooring_data.CO2.DateTime(i))<1
                idx_lat = knnsearch(getfield(ERA5,fields{s},'Lat'),mooring_data.CO2.Latitude(i));
                idx_lon = knnsearch(getfield(ERA5,fields{s},'Lon'),mooring_data.CO2.Longitude(i));

                mooring_data.CO2.ERA5_wsp(i)=sqrt((getfield(ERA5,fields{s},'u_10',{idx_lon,idx_lat,idx}))^2 + (getfield(ERA5,fields{s},'v_10',{idx_lon,idx_lat,idx}))^2);
                mooring_data.CO2.ERA5_time(i)=getfield(ERA5,fields{s},'Time',{idx});
                mooring_data.CO2.ERA5_lat(i)=getfield(ERA5,fields{s},'Lat',{idx_lat});
                mooring_data.CO2.ERA5_lon(i)=getfield(ERA5,fields{s},'Lon',{idx_lon});
                mooring_data.CO2.ERA5_airtemp(i)=getfield(ERA5,fields{s},'t_2',{idx_lon,idx_lat,idx});
            end
        end
    end
end


%%%% now calculate fluxes with pCO2
%%% result in mmol m-2 d-1
for i=1:height(mooring_data.CO2)
    if ~isnan(mooring_data.CO2.dpCO2(i))
        if mooring_data.CO2.ERA5_lat(i) == 0
            mooring_data.CO2.flux_mmol_m2_d1(i)=FCO2_CWE(mooring_data.CO2.dpCO2(i),mooring_data.CO2.SST_C_(i),mooring_data.CO2.Salinity(i),mooring_data.CO2.wsp(i));
        else
            mooring_data.CO2.flux_mmol_m2_d1(i)=FCO2_CWE(mooring_data.CO2.dpCO2(i),mooring_data.CO2.SST_C_(i),mooring_data.CO2.Salinity(i),mooring_data.CO2.ERA5_wsp(i));
        end
    else
        mooring_data.CO2.flux_mmol_m2_d1(i)=NaN;
    end
end

mooring_data.CO2.flux_mol_m2_y1 = (mooring_data.CO2.flux_mmol_m2_d1/1000)*365;


clearvars -except mooring_data raw_data

path =('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts\Argo-air-sea-flux-paper');
cd(path)
save('complete_SOFS_data.mat')

% visualize the data

figure()
yyaxis left
plot(mooring_data.CO2.DateTime, mooring_data.CO2.dpCO2, 'ob')
ylabel('dpCO2 - uatm')
yyaxis right
plot(mooring_data.CO2.DateTime, mooring_data.CO2.flux_mol_m2_y1, '*r')
ylabel('air sea flux - mol m^-^2 y^-^1')
xlabel('Time')
legend('dpCO_2' , 'air sea flux')
hold on
yline(0,'HandleVisibility','off')


figure()
plot(mooring_data.CO2.wsp_time, mooring_data.CO2.wsp,'or')
hold on
plot(mooring_data.CO2.ERA5_time, mooring_data.CO2.ERA5_wsp,'ob')
hold off
xlabel('Time')
ylabel('wsp m/s')
legend('mooring' , 'ERA5')

figure()
hold on
plot(mooring_data.CO2.DateTime(mooring_data.CO2.Salinity~=-999), mooring_data.CO2.Salinity(mooring_data.CO2.Salinity~=-999), '+b')
plot(mooring_data.CO2.DateTime(mooring_data.CO2.SST_C_~=-999), mooring_data.CO2.SST_C_(mooring_data.CO2.SST_C_~=-999), '+k')
hold off
xlabel('Time')
legend('PSAL', 'SST')

