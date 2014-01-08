function [M,R] = CPRBayes_dice(data,distr,SETTINGS)
%CPRBAYES_DICE Summary of this function goes here
%   Detailed explanation goes here

x = (1:SETTINGS.LEN)';
seg = floor((x-1)./(SETTINGS.LEN./(SETTINGS.DICE)))+1;
M = SETTINGS.PRIOR_M;
R = zeros(length(M),1);
for i = 1:SETTINGS.DICE
    subdex = find(seg==i);
        subargin = [SETTINGS.VARARGIN 'timestamps' SETTINGS.TIME(subdex) 'messages' 'off' 'pruning' 0 'dice' 1 'depth' 1];
    [Md,~,statd] = CPRBayes(data(subdex,:),distr,subargin{1,:});
    if statd.post_ratios(2) > SETTINGS.THRESH;
        M = [M;Md(2)+subdex(1)];
        R = [R;statd.post_ratios(2)];
    end
end
M = sort(M);

end

