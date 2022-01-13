load('SOTS_float_data.mat')

%%%%% set the reference depth (1500m as per Williams, 2017)
ref_depth = ones(1,size(SOTS_float_data.pH,2))*1500;
%%% set the tolerance to the reference depth
ref_depth_t = 100;
%%% find the index of each float profile where the depth is closest to the
%%% reference depth
for i=1:size(SOTS_float_data.pH,2)
    pres_ref_id(i) = knnsearch(SOTS_float_data.pres(:,i), ref_depth(i));
    pres_ref(i) = SOTS_float_data.pres(pres_ref_id(i),i);
end

%%% take out the depth that are not within the depth tolerance
pres_ref(pres_ref >= (ref_depth(1) + ref_depth_t))=NaN;
pres_ref(pres_ref <= (ref_depth(1) - ref_depth_t))=NaN;

%%% find the profile numbers where we have a deep enough measurement
pH_ref_profile_ind=find(~isnan(pres_ref));

pH_ref_depth_ind=pres_ref_id(pH_ref_profile_ind);

%%% find the pH at the reference depths
for i=1:length(pH_ref_profile_ind)
    pH_ref.pH(i) = SOTS_float_data.pH(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.temp(i) = SOTS_float_data.TEMP(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.pres(i) = SOTS_float_data.pres(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.psal(i) = SOTS_float_data.psal(pH_ref_depth_ind(i),pH_ref_profile_ind(i));
    pH_ref.profile_no(i) = pH_ref_profile_ind(i);
    %%%% calculate the pH at 25C and 0dbar
    [pH_25] = CO2SYS(2290,pH_ref.pH(i),1,3,pH_ref.psal(i),...
        pH_ref.temp(i),25, pH_ref.pres(i),0,2.8,0.9,2,0,1,10,1,2,2);
    pH_ref.pH_25(i) = pH_25(21);
    %%% calculate the bias correction
    pH_ref.pH_corr(i) = -0.034529 * pH_ref.pH_25(i) + 0.26709;
end

%%%% now each pH profile has to be corrected by adding the bias
%%% I will use the correction that is closest to each profile
for i = 1:size(SOTS_float_data.pH,2)
%for i=15   
    ind = knnsearch(pH_ref.profile_no(:),i);
    pH_LD_corr(:,i) = SOTS_float_data.pH_LIR_Deep(:,i)+pH_ref.pH_corr(ind);
    pH_LS_corr(:,i) = SOTS_float_data.pH_LIR_Shallow(:,i)+pH_ref.pH_corr(ind);
    pH_LD_corr(:,i) = SOTS_float_data.pH_Williams_Deep(:,i)+pH_ref.pH_corr(ind);
end



