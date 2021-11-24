% this script calculates SOFS mooring air-sea fluxes

% download data from
% https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0118546/
% stitch CSV files together,
% calculate flux from dfCO2, SST and S in the files
% find windspeed on thredd server (from the mooring)
% translate windspeed from whichever mooring height to 10m
% with Sutton, 2017 equation
% u10[m s-1] = u [m s-1] / (1 + (sqrt(0.0011)/0.4) * ln(Z/10))
% Z is wind anemometer height [m]