% this script takes the monthly Cape Grim xCO2 (ppm) in dry air
% measurements and converts them into pCO2 in air (with the actual
% water vapour pressure calculated). The data is then interpolated
% with a piecewise cubic Hermite spline (as per Gray paper)
file_path = 'C:\Users\cawynn\cloudstor\Air sea flux manuscript\';

file = [file_path 'CapeGrim_CO2_data_download_08112021.xlsx'];

function [CGdata, CGdata_interp] = Cape_Grim_calc(file)

    CGdata = readtable(file, 'Sheet', 'CapeGrim_CO2_data_download_0811');
    % clean out the NAN lines (which were text in the original)
    CGdata(isnan(CGdata.YYYY),:)=[];

    CGdata.datestr = strcat(num2str(CGdata.YYYY),{'-'},num2str(CGdata.MM),{'-'}, num2str(CGdata.DD));

    date = datetime(CGdata.datestr, 'InputFormat','yyyy-MM-dd');

    CGdata_interp = [];
    CGdata_interp.date = date;
    CGdata_interp.ppm = CGdata.CO2_ppm_;
    CGdata_interp.datehours=(datenum(date)-datenum('1990-1-1'))*24;
    xq2 = CGdata_interp.datehours(1):1:CGdata_interp.datehours(end);
    CGdata_interp.ppm_spl = pchip(CGdata_interp.datehours,CGdata.CO2_ppm_,xq2);
    %p = pchip(CGdata_interp.datehours,CGdata.CO2_ppm_,xq2);
end



% figure()
% scatter(CGdata_interp.datehours,CGdata.CO2_ppm_,'o')
% hold on
% scatter(xq2,CGdata_interp.ppm_spl,'*')
% hold off



