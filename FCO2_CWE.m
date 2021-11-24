

function [F_CO2]=FCO2(dpCO2,T,S,u)
   
    %%Function for the calculation or air-sea CO2 flux
        %   modified after Cecilia Chapa Balcorta 
        %   Ensenada, Baja California
        %   UNIVERSIDAD AUTÓNOMA DE BAJA CALIFORNIA/ UNIVERSIDAD DEL MAR
    
        %   filename: FCO2.m
        %   created: 07/04/2014
        %   modified: 03/23/2015
        %   CWE modified: 23/11/2021
        %   copyright @ 2015 Cecilia Chapa-Balcorta
        
       
        %INPUT
    
   %pCO2_agua= seawater pCO2 (uatm)
   %pCO2_atm=  atmospheric pCO2 (uatm)
   %T=  Temperature (Celsius)
   %S=  Salinity 
   %u = Wind speed (m/s)
   
    
   %Air-sea CO2 is calculated as follows:
   
           % FCO2 =K*a(dpCO2) 
    
    %Where
    %K=is the transfer velocity according to Wanninkhof (2014).
    %a = CO2 solubility constant according to Weiss (1974)
    %dpCO2 is the difference of air and seawater pCO2 
    
    
    %%%%% CO2 Transfer velocity calculation %%%%%%%%%
       
   Sc=Schmidt(T);
%    if u<=6;    K=0.31*(u.^2).*((Sc./660).^-0.5); %for slower steadier wind
%    %K=K*3600*24;  %conversion to   m/day
%    else    K=0.39*(u.^2).*((Sc./660).^-0.5);
%    end
   % k units in cm h-1, u units m s-1
   K = 0.251 *(u.^2).*((Sc./660).^-0.5); % updated as per Wanninkhof, 2014
   % 6-hourly windproduct at 0.25 degree resolution, appropriate for 3-15 m s-1
   
   %dpCO2=pCO2_agua-pCO2_atm;     %%calculation of delta pCO2
   a=Ko_weiss(T,S);    %Solubility in mmol L^-1 atm^-1 or mmol m^-3 uatm^-1
    
F_CO2 =0.24*K.*a.*dpCO2; %CO2 flux (mmol m^-2 d^-1)
end

%%%%%Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%****Solubility constant (Weiss, 1974) ***********************************

function [Ko]=Ko_weiss(T,S)
%
A=[-60.2409, 93.4517, 23.3585];  %mol/Kg.atm
B=[0.023517, -0.023656, 0.0047036]; %mol/Kg.atm
T=T+273.15; %Conversion from Celsius degrees to Kelvins
Ln_Ko=A(1)+(A(2).*(100./T))+(A(3).*log(T./100))+S.*(B(1)+(B(2).*(T./100))+(B(3).*(T./100).^2));
Ko=exp(Ln_Ko);

end

%******** Schmidt Number*********
% updated as per Wanninkhof, 2014
    %For water of salinity=35 and temperature range 0-30°C    %%%%%%%%%%%%%
    
function [Sc]=Schmidt(T)
        A = 2116.8;     B = -136.25;     C = 4.7353;     D = -0.092307;
        E = 0.0007555;
        Sc= A + (B.*T)+(C.*T.^2)-(D.*T.^3)+(E.*T.^4);
end
        
