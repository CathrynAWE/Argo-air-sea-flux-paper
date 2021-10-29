% SOTS 2020
files{1} = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V09 SOTS/IMOS_SOOP-CO2_GST_20200827T030011Z_VLMJ_FV01.nc';
% SOLACE
files{2} = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2020_V08 SOLACE/IMOS_SOOP-CO2_GST_20201204T043818Z_VLMJ_FV01.nc';
% TEMPO
files{3} = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2021_V01 TEMPO/IMOS_SOOP-CO2_GST_20210129T004303Z_VLMJ_FV01.nc';
% SOTS 2021
files{4} = 'C:/Users/cawynn/cloudstor/Air sea flux manuscript/IN2021_V02 SOTS/IMOS_SOOP-CO2_GST_20210414T121956Z_VLMJ_FV01.nc';

fig = figure()
for i = 1:length(files)
    pCO2_sw = ncread(files{i},'fCO2SW_UATM');
    pCO2_atm = ncread(files{i},'fCO2ATM_UATM_INTERPOLATED');
    T = ncread(files{i},'TEMP');
    S = ncread(files{i},'PSAL');
    u = ncread(files{i},'WSPD');
    time = ncread(files{i}, 'TIME') + datetime(1950,1,1);
    lat = ncread(files{i}, 'LATITUDE');
    lon = ncread(files{i}, 'LONGITUDE');

    [F_CO2, dpCO2]=FCO2(pCO2_sw, pCO2_atm,T,S,u);

    
    s=5;
    scatter3(time(F_CO2>-100),lat(F_CO2>-100),lon(F_CO2>-100),[],F_CO2(F_CO2>-100),'filled')
    c = colorbar;
    c.Label.String = 'Oxygen';
    hold on
    xlabel('Time')
    ylabel('Latitude')
    zlabel('Longitude')
end
