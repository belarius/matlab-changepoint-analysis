function [M,P,stats] = CPRBayesForward(D,distr,varargin)
%CPRBAYESFORWARD Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    error('CPRBayes:TooFewInputs','Too few inputs; Both the observations and the distribution must be specified');
end
distr = lower(distr);
if ismember(distr,{'binomial','geometric','poisson','exponential','linear','normal','uniform','multiple linear','multivariate normal','multinomial'})==0
    error('CPRBayes:UnknownDistr','Distribution not recognized');
end

if ~isempty(varargin)
    if iscell(varargin{1})
        varargin = varargin{1};
    end
end

varargin = ['thresh' 20 varargin];

data = SetupData(D,distr);
SETTINGS = SetupSettings(data,distr,varargin);

len = SETTINGS.LEN;
ops = SETTINGS.OPS;
p_c = SETTINGS.PRIOR_C;
M = SETTINGS.PRIOR_M;
times = SETTINGS.TIME;
a = SETTINGS.ALPHA;
regressors = SETTINGS.REGRESSORS;
thresh = SETTINGS.THRESH;
weight = SETTINGS.WEIGHT;

if SETTINGS.CONF
    Conf_Placeholder = SETTINGS.CONF;
    SETTINGS.CONF = 0;
end

R = nan(length(M),1);
k = nan(size(D,1),1);
td = nan(size(D,1),1);

ticID = tic;

mdex = 1;
for i = 2:len
    subargin = [varargin 'messages' 'off' 'priorc' p_c(1)./p_c(2) 'alpha' {a} 'timestamps' times(M(mdex)+1:i)];
    [Mp,~,statp] = CPRBayes(data(M(mdex)+1:i,:),distr,subargin);
    if length(Mp) > 2
        sortblock = sortrows([M R;Mp(2)+M(mdex) statp.post_ratios(2)],1);
        M = sortblock(:,1);
        R = sortblock(:,2);
        k(M(mdex)+1:i) = [statp.seg_numer(:,1);nan];
        td(M(mdex)+1:i) = [statp.seg_denom(:,1);nan];
        if length(M) > 3
            subargin = [varargin 'messages' 'off' 'priorc' p_c(1)./p_c(2) 'alpha' {a} 'timestamps' times(M(mdex-1)+1:M(mdex+1))];
            [Mp,~,statp] = CPRBayes(data(M(mdex-1)+1:M(mdex+1),:),distr,subargin);
            if length(Mp) == 2
                M(mdex) = [];
                R(mdex) = [];
                k(M(mdex-1)+1:M(mdex)) = [statp.seg_numer(:,1);nan];
                td(M(mdex-1)+1:M(mdex)) = [statp.seg_denom(:,1);nan];
            else
                M(mdex) = Mp(2)+M(mdex-1);
                R(mdex) = statp.post_ratios(2);
                k(M(mdex-1)+1:M(mdex)) = [statp.seg_numer(1:Mp(2)-1,1);nan];
                td(M(mdex-1)+1:M(mdex)) = [statp.seg_denom(1:Mp(2)-1,1);nan];
                mdex = mdex+1;
            end
        else
            mdex = mdex+1;
        end
        p_c(1) = length(M)-1;
    end
end
k(M(mdex)+1:end) = [statp.seg_numer(:,1);nan];
td(M(mdex)+1:end) = [statp.seg_denom(:,1);nan];
k(end) = []; td(end) = [];

% mdex = 1;
% subargin = [varargin 'messages' 'off' 'priorc' p_c(1)./p_c(2) 'alpha' a 'depth' 1];
% while length(M) > mdex+1
%     [Mp,~,statp] = CPRBayes(data(M(mdex)+1:M(mdex+2),:),distr,subargin);
%     if statp.post_ratios < statp.threshold
%         M(mdex+1) = [];
%     else
%         M(mdex+1) = M(mdex)+Mp(2);
%         mdex = mdex+1;
%     end
% end

%==CONFIDENCE INTERVALS
SETTINGS.CONF = Conf_Placeholder;
if SETTINGS.CONF
    C = EstimateCI(data,distr,M,SETTINGS);
else
    C = [];
end

%===ESTIMATE PARAMETERS===
P = EstimateParameters(data,distr,M,SETTINGS.KNOWN_P,SETTINGS.ALPHA,SETTINGS.PARAMS);

%===CLOSING DOWN===
ElapsedTime = toc(ticID);
if SETTINGS.MESSAGES == 1
    disp(['Elapsed time is ' num2str(ElapsedTime) ' seconds.'])
end

stats = struct('model',{M},'distr',{distr},'params',{P},'timestamps',times,'regressors',regressors,'weights',weight,'threshold',{thresh},'post_ratios',{R},'seg_numer',k,'seg_denom',td,'conf',SETTINGS.CONF_A,'conf_int',{C},'prior_change',{p_c(1)./p_c(2)},'prior_model',{SETTINGS.ALPHA},'elapsed_time',ElapsedTime);

end