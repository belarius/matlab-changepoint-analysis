function SETTINGS = SetupSettings(data,distr,varargin)
%SETUPSETTINGS Summary of this function goes here
%   Detailed explanation goes here

len = size(data,1);
ops = size(data,2);
t = (1:len)';

nanv = find(sum(isnan(data),2)==1);
if ~isempty(nanv)
    warning('CPRBayes:NaNDetected','NaN values detected at the following indices: %s',nanv);
    t(nanv) = [];
    data(nanv,:) = [];
    len = size(data,1);
end

M = [0;len];
SETTINGS = struct('LEN',len','OPS',ops,'PRIOR_M',M,'ALPHA',[],'DEPTH',0,'KNOWN_P',[],'PRIOR_C',[1 len-1],'THRESH',10,'TIME',t,'REGRESSORS',[ones(len,1) t],'DICE',1,'WEIGHT',ones(len-1,1),'MESSAGES',1,'PRUNING',0,'CONF',1,'CONF_A',0.95,'PARAMS',0,'VARARGIN',varargin);

if ~isempty(varargin)
    if iscell(varargin{1})
        varargin = varargin{1};
    end
end

%===PARAMETER INTERPRETATION===
if mod(length(varargin),2)
    error('CPRBayes:BadlyFormedParameters','Incorrectly formatted parameters');
else
    rsargin = reshape(varargin,2,length(varargin)./2)';
    for i = 1:length(rsargin(:,1))
        switch rsargin{i,1}
            case 'alpha'
                SETTINGS.ALPHA = rsargin{i,2};
            case 'conf'
                SETTINGS.CONF_A = rsargin{i,2};
            case 'depth'
                SETTINGS.DEPTH = rsargin{i,2};
            case 'knownp'
                SETTINGS.KNOWN_P = rsargin{i,2};
                if sum(isnan(SETTINGS.KNOWN_P)) == length(SETTINGS.KNOWN_P)
                    error('CPRBayes:NoFreeParametersRemaining','You must have at least one unknown parameter, designated as NaN');
                end
            case 'knownM'
                SETTINGS.PRIOR_M = rsargin{i,2};
            case 'priorc'
                SETTINGS.PRIOR_C = [rsargin{i,2}.*len len];
            case 'timestamps'
                trash = rsargin{i,2};
                trash(nanv) = [];
                SETTINGS.TIME = trash;
                SETTINGS.REGRESSORS(:,2) = SETTINGS.TIME;
            case 'regressors'
                trash = rsargin{i,2};
                trash(nanv,:) = [];
                SETTINGS.REGRESSORS = trash;
            case 'thresh'
                SETTINGS.THRESH = rsargin{i,2};
            case 'dice'
                SETTINGS.DICE = rsargin{i,2};
            case 'weight'
                SETTINGS.WEIGHT = rsargin{i,2};
            case 'messages'
                switch rsargin{i,2}
                    case 'on'
                        SETTINGS.MESSAGES = 1;
                    otherwise
                        SETTINGS.MESSAGES = 0;
                end
            case 'pruning'
                SETTINGS.PRUNING = rsargin{i,2};
            case 'conftype'
                switch rsargin{i,2}
                    case 'itemwise'
                        SETTINGS.CONFTYPE = 1;
                    case 'none'
                        SETTINGS.CONFTYPE = 0;
                end
            case 'paramtype'
                switch rsargin{i,2}
                    case 'bayesian'
                        SETTINGS.PARAMS = 0;
                    case 'frequentist'
                        SETTINGS.PARAMS = 1;
                    case 'robust'
                        SETTINGS.PARAMS = 2;
                end
        end
    end
   for i = 1:length(rsargin(:,1))
        switch rsargin{i,1}
            case 'regressors'
                SETTINGS.REGRESSORS = rsargin{i,2};
        end
   end
end

if isempty(SETTINGS.KNOWN_P)
    SETTINGS.KNOWN_P = DefaultParams(distr,ops);
end

if isempty(SETTINGS.ALPHA)     %Each distribution has its own default alpha array, corresponding to a weak assumption with a proper integral
    SETTINGS.ALPHA = DefaultAlpha(data,distr,SETTINGS.KNOWN_P,ops,SETTINGS.REGRESSORS);
end

if ~strcmp(distr,'linear')&&~strcmp(distr,'multiple linear')
    SETTINGS.REGRESSORS = NaN;
end

if SETTINGS.DICE < 1
    SETTINGS.DICE = 1;
else
    SETTINGS.DICE = floor(SETTINGS.DICE);
end


end

