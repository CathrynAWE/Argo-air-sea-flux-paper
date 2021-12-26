% get ftp data from ifremer ARGO server
% CS, 3.6.2019

% connect to server
ftpobj = ftp('ftp.ifremer.fr');
%ftpobj = ftp('ftp.usgodae.org');

% go to the folder with all the float files
% cd(ftpobj,'ifremer/argo/dac/csiro/5906623'); % SOTS
cd(ftpobj,'ifremer/argo/dac/csiro/5906624'); % 55S


% find the BGC files or find the Sprof
%BGC = dir(ftpobj,'B*.*');
BGC = dir(ftpobj,'*Sprof.nc');

% cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\5906623')
cd('C:\Users\cawynn\cloudstor\Air sea flux manuscript\BGC_Argo\5906624')


% and loop through them to get out the data
for k=1:length(BGC)
    
    nprof = BGC(k).name;
    
    mget(ftpobj,nprof);
  
    
end
