%% convert from mole fraction and do CO2sys computations on SOTS data

% load the file from Pete
%file='pCO2-merge.nc';
file='pCO2-merge-2021-08-15.nc';

time = ncread(file, 'TIME') + datetime(1950,1,1);
psal = ncread(file, 'PSAL');
temp = ncread(file, 'TEMP');

xCO2_SW = ncread(file, 'xCO2_SW');
xCO2_AIR = ncread(file, 'xCO2_AIR');

pressure = ncread(file, 'pressure');

d = time;
doy = day(d,'dayofyear');

sss = psal;
sst = temp;

% convert from mole fraction to pCO2
slp = pressure./101.325;
svp = exp(25.4543-67.4509.*(100./(sst+273.15)) - 4.8489.*log((sst+273.15)./100) - 0.000544.*sss);
pco2 = xCO2_SW.*slp.*(1-svp);
pco2atm = xCO2_AIR.*slp.*(1-svp);
Dpco2 = pco2 - pco2atm;

%% do the CO2 system computations
alk = 39.23.*sss + 937.3;

par1type =    4; % The first parameter supplied is of type "4", which is "pCO2"
par1     = pco2; % value of the first parameter
par2type =    1; % The second parameter supplied is of type "1", which is "alkalinity"
par2     = alk; % value of the second parameter, 
sal      =   sss; % Salinity of the sample
tempin   =   sst; % Temperature at input conditions
presin   =   1; % Pressure    at input conditions
tempout  =    0; % Temperature at output conditions - doesn't matter in this example
presout  =    0; % Pressure    at output conditions - doesn't matter in this example
sil      =   4; % Concentration of silicate  in the sample (in umol/kg)
po4      =    1; % Concentration of phosphate in the sample (in umol/kg)
pHscale  =   1; % pH scale at which the input pH is reported ("1" means "Total Scale")  - doesn't matter in this example
k1k2c    =    10; % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("10" means "Lueker")
kso4c    =    1; % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
sofs_co2sys=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

dic = sofs_co2sys(:,2);
pH = sofs_co2sys(:,3);
omega = sofs_co2sys(:,16);

clear par1 par2 par1type par2type sal tempin tempout presin presout sil po4 pHscale k1k2c kso4c

