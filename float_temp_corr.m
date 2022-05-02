%%%%% set the reference depth (1500m as per Williams, 2017)
%%%%% I will use this for the LD and WD corrected profiles
ref_depth_D = 1500;
%%% set the tolerance to the reference depth
ref_depth_t_D = 50;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i= 1:size(SOTS_float_data.pH,2)
    pres_ref_id_D = SOTS_float_data.pres(:,i)>=(ref_depth_D-ref_depth_t_D) & SOTS_float_data.pres(:,i)<=(ref_depth_D+ref_depth_t_D); 
    %%% find the TEMP at the reference depths
    pH_ref.temp_D(1:length(SOTS_float_data.TEMP),i) = mean(SOTS_float_data.TEMP(pres_ref_id_D,i),'omitnan');
 
end

%%% find the profile numbers where we have a deep enough measurement
pH_ref.profile_no_D = find(~isnan(pH_ref.temp_D(1,:)));

       
%%%% now the temp correction for offset and drift is calculated for each
%%%% profile
%%% I will use the temp_ref that is closest to each profile
for i = 1:size(SOTS_float_data.pH,2)
% for i=6   
    ind = knnsearch(pH_ref.profile_no_D(:),i);
    idx = pH_ref.profile_no_D(ind);
    % calculate the temp corrections, Temp needs conversion from C to K
    pH_ref.temp_corr_D(:,i) = (pH_ref.temp_D(:,idx)+273.16)./(SOTS_float_data.TEMP(:,i)+273.16);
end

