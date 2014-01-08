function C = EstimateCI(data,distr,M,SETTINGS)
%ESTIMATECI Summary of this function goes here
%   Detailed explanation goes here

times = SETTINGS.TIME;
weight = SETTINGS.WEIGHT;

C = nan(length(M),2);
if length(M) > 2
    for m = 2:length(M)-1
        SUBSETTINGS = SETTINGS;
        SUBSETTINGS.LEN = M(m+1)-M(m-1);
        SUBSETTINGS.PRIOR_C(2) = SUBSETTINGS.LEN;
        SUBSETTINGS.PRIOR_M(2) = SUBSETTINGS.LEN;
        SUBSETTINGS.TIME = times(M(m-1)+1:M(m+1));
        SUBSETTINGS.WEIGHT = weight(M(m-1)+1:M(m+1)-1);
        if strcmp(distr,'linear')
            SUBSETTINGS.ALPHA{5} = SETTINGS.ALPHA{5}(M(m-1)+1:M(m+1),:);
        end
        subd = data(M(m-1)+1:M(m+1),:);
        [~,ci,~,~] = CPRBayes_single(subd,distr,SUBSETTINGS,SUBSETTINGS.TIME,SUBSETTINGS.WEIGHT);
        C(m,:) = ci;
    end
end

end

