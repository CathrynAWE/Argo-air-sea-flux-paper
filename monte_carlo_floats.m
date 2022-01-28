% this script calculates the monte carlo simulations for air sea flux
% calculations from the floats.
clear all

addpath('C:\Users\cawynn\cloudstor\Air sea flux manuscript\Matlab scripts');


% float_ID = '5906623' % SOTS float
float_ID = '5906624' % 55S float

% change structure name!
monte_carlo_S55 = [];

% number of simulations
no_runs = 1000;

if float_ID == '5906623' % SOTS float
    load('SOTS_float_data.mat')
    data = SOTS_float_data;
elseif float_ID == '5906624' % 55S float
    load('S55_float_data.mat')
    data = S55_float_data;
end

% create the monte carlo structure, load the necessary variables into it

monte_carlo_S55.raw.pCO2_LD_corr = data.pCO2_L_D_corr;
monte_carlo_S55.raw.pCO2_LS_corr = data.pCO2_L_S_corr;
monte_carlo_S55.raw.pCO2_uatm = data.pCO2_uatm;
monte_carlo_S55.raw.swtemp = data.TEMP_20;
monte_carlo_S55.raw.psal = data.PSAL_20;
monte_carlo_S55.raw.wsp = data.wsp;
monte_carlo_S55.raw.time = data.time;
monte_carlo_S55.raw.lat = data.lat;
monte_carlo_S55.raw.lon = data.lon;
monte_carlo_S55.raw.month = data.Month;
monte_carlo_S55.raw.year = data.year;
monte_carlo_S55.raw.flux_LD_corr = data.flux_L_D_corr;
monte_carlo_S55.raw.flux_LS_corr = data.flux_L_S_corr;
monte_carlo_S55.raw.flux_LD_moav = data.flux_LD_corr_mo_ave;
monte_carlo_S55.raw.flux_LS_moav = data.flux_LS_corr_mo_ave;


% define the uncertainties 
% randomly assigned --> pCO2 * uncertainty + pCO2
%randi([0, 1], 1) * normrnd(0,0.022);
% systematically assigned to all --> pCO2 * uncertainty + pCO2
%normrnd(0,0.018);

% atmospheric pCO2 --> pCO2 + uncertainty
%normrnd(0,0.15);

% temperature --> temp + uncertainty
pd1 = makedist('Uniform','lower',-0.002,'upper',0.002);
%random(pd1);

% salinity --> psal + uncertainty
pd2 = makedist('Uniform','lower',-0.01,'upper',0.01);
%random(pd2);

% windspeed --> wsp + wps * uncertainty
%normrnd(0,1.5);

% now add the uncertainties
for i = 1:size(data.time,2)
    for n = 1:no_runs
        
        monte_carlo_S55.simulation.pCO2_LD_corr(n,i) = ...
            data.pCO2_L_D_corr(i) + (data.pCO2_L_D_corr(i)*(randi([0, 1], 1) * normrnd(0,0.022)))+ ...
            (data.pCO2_L_D_corr(i)*(normrnd(0,0.018)));

        monte_carlo_S55.simulation.pCO2_LS_corr(n,i) = ...
            data.pCO2_L_S_corr(i) + (data.pCO2_L_S_corr(i)*normrnd(0,0.018))+ ...
            (data.pCO2_L_S_corr(i)*normrnd(0,0.018));

        monte_carlo_S55.simulation.pCO2_uatm(n,i) = data.pCO2_uatm(i) + normrnd(0,0.15);

        monte_carlo_S55.simulation.swtemp(n,i) = data.TEMP_20(i) + random(pd1);
        monte_carlo_S55.simulation.psal(n,i) = data.PSAL_20(i) + random(pd2);
        monte_carlo_S55.simulation.wsp(n,i) = data.wsp(i) + normrnd(0,1.5);
    end
end


for i = 1:size(data.time,2)
    for n = 1:no_runs
    
        [FCO2_LD] = FCO2_CWE_monte_carlo(monte_carlo_S55.simulation.pCO2_LD_corr(n,i),...
            monte_carlo_S55.simulation.pCO2_uatm(n,i),monte_carlo_S55.simulation.swtemp(n,i),...
            monte_carlo_S55.simulation.psal(n,i),monte_carlo_S55.simulation.wsp(n,i));
        monte_carlo_S55.simulation.flux_LD(n,i) = [FCO2_LD];
        
        [FCO2_LS] = FCO2_CWE_monte_carlo(monte_carlo_S55.simulation.pCO2_LS_corr(n,i),...
            monte_carlo_S55.simulation.pCO2_uatm(n,i),monte_carlo_S55.simulation.swtemp(n,i),...
            monte_carlo_S55.simulation.psal(n,i),monte_carlo_S55.simulation.wsp(n,i));
        monte_carlo_S55.simulation.flux_LS(n,i) = [FCO2_LS];
        
    end
end

% I will only use profiles that happened in 2021

% subscript for accumarray
n = find(data.year ==2021,1);
for i = 1:12
    c=find(data.Month(n:size(data.time,2))==i);
    c1=c(1);
    c2=c(end);
    monte_carlo_S55.simulation.flux_LD_moav(i)=...
        mean(monte_carlo_S55.simulation.flux_LD(:,c1:c2),'all','omitnan');
    a=mat2cell(monte_carlo_S55.simulation.flux_LD(:,c1:c2),no_runs,length(c1:c2));
    monte_carlo_S55.simulation.flux_LD_moSD(i)=...
        cellfun(@(x) std(x,0,[1,2],'omitnan'),a,'UniformOutput',true);
    
    monte_carlo_S55.simulation.flux_LS_moav(i)=...
        mean(monte_carlo_S55.simulation.flux_LS(:,c1:c2),'all','omitnan');
    a=mat2cell(monte_carlo_S55.simulation.flux_LS(:,c1:c2),no_runs,length(c1:c2));
    monte_carlo_S55.simulation.flux_LS_moSD(i)=...
        cellfun(@(x) std(x,0,[1,2],'omitnan'),a,'UniformOutput',true);
end


if float_ID == '5906623'
    clearvars -except monte_carlo_SOTS
    save('monte_carlo_SOTS.mat')
elseif float_ID == '5906624'
    clearvars -except monte_carlo_S55
    save('monte_carlo_S55.mat')
end

