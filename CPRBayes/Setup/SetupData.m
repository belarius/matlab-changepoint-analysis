function data = SetupData(D,distr)
%SETUPDATA Summary of this function goes here
%   Detailed explanation goes here

switch distr
    case {'binomial','multinomial'}
        ops = size(D,2);
        if sum(~(rem(D(:),1) == 0)) > 0
            error('CPRBayes:NotDiscrete','data() must be discrete if using the %s model',distr);
        end
        if strcmp('binomial',distr)
            if ops == 1
                dex = find(~((D==1)|(D==0)));
                if ~isempty(dex)
                    error('CPRBayes:NotBernoulli','When data() is binomial and has one column, it must consist of only ones and zeros');
                end
                data = [D 1-D];
            elseif ops == 2
                data = [D(:,1) D(:,2)-D(:,1)];
                if ~isempty(find(data(:,2)<0))
                    error('CPRBayes:NotSuccessVsTrial','When data() is binomial and has two columns, data(:,2) must >= data(:,1)');
                end
            else
                error('CPRBayes:UniVsMultivariate','data() has more than two outcomes; consider using the multinomial distribution instead');
            end
        elseif ops == 1 && strcmp('multinomial',distr)
            error('CPRBayes:UniVsMultivariate','data() has only 1 column; consider using the binomial distribution instead');
        else
            data = D;
        end
    case {'geometric','poisson'}
        if sum(~(rem(D(:),1) == 0)) > 0
            error('CPRBayes:NotDiscrete','data() must be discrete if using the %s model',distr);
        end
        if size(D,2) > 1
            error('CPRBayes:UniVsMultivariate','data() may only have 1 column; multivariate %s not supported',distr);
        end
        if sum(mod(D,1)) > 0
            warning('CPRBayes:DataDistributionSupportMismatch','%s data may only contain integer values; rounding data() to the nearest integer',distr);
        end
        data = round(D);
    case {'exponential','linear','normal','uniform'}
        if size(D,2) > 1
            error('CPRBayes:UniVsMultivariate','data() may only have 1 column; multivariate %s not supported',distr);
        end
        data = D;
    case {'multiple linear','multivariate normal'}
        if size(D,2) == 1
            error('CPRBayes:UniVsMultivariate','data() has only 1 column; consider using a univariate distribution instead');
        end
        data = D;
end
end
