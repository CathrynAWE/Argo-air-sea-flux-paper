% get ftp data from ifremer ARGO server
% CS, 3.6.2019

%float_ID = '5904686'; %55S complimentary float
%  float_ID = '5906623'; % SOTS
float_ID = '5906624'; % 55S


% connect to server
ftpobj = ftp('ftp.ifremer.fr');
%ftpobj = ftp('ftp.usgodae.org');

% go to the folder with all the float files
float_ifremer_path = ['ifremer/argo/dac/csiro/' float_ID];
cd(ftpobj, float_ifremer_path);
% cd(ftpobj,'ifremer/argo/dac/csiro/5906623'); % SOTS
% cd(ftpobj,'ifremer/argo/dac/csiro/5906624'); % 55S
% cd(ftpobj,'ifremer/argo/dac/csiro/5904686'); %55S complimentary float

% find the BGC files or find the Sprof
%BGC = dir(ftpobj,'B*.*');
BGC = dir(ftpobj,'*Sprof.nc');

float_storage_path = ['C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\' float_ID];
cd(float_storage_path);
% cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\5906623')
% cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\5906624')
% cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\5904686'); %55S complimentary float

% and loop through them to get out the data
for k=1:length(BGC)
    
    nprof = BGC(k).name;
    
    mget(ftpobj,nprof);
  
    
end
