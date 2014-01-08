function a = DefaultAlpha(data,distr,knownp,ops,regressors)
%DEFAULTALPHA Sets the prior hyperparameters according to the default empirical Bayes 'rule of thumb.'
%   ~fill in later~
switch distr
    case {'binomial','geometric'}
        a = [0.5 0.5];
    case 'exponential'
        a = [1 median(data)];
    case 'poisson'
        a = [mean(data) 1];
    case 'multinomial'
        a = ones(1,ops).*0.5;
    case 'linear'
        len = length(data(:,1));
        a = cell(1,5);
        a{5} = regressors;
        a{1} = regress(data,a{5});
        a{2} = 1;
        a{3} = 1;
        a{4} = (a{5}'*a{5})./(len.^2);
    case 'normal'
        if isnan(knownp(1)) && isnan(knownp(2))
            a = [median(data) 0.1 0.1 mad(data,1)*1.4826]; %Median(x) and MedianAbsoluteDeviation(x) used to estimate mean and SD
        elseif isnan(knownp(1))
            a = [median(data) 0.1]; %Median(x) used to estimate mean
        else
            a = [0.1 mad(data,1)*1.4826]; %MedianAbsoluteDeviation(x) used to estimate SD
        end
    case 'uniform'  %The fourth alpha is a gap delimiter that depends on a(2) and a(3) but can be defined independently.
        df = diff(sort(data));
        df(df==0) = [];
        if isempty(df)
            df = 1;
        end
        a = [1 knownp(1) knownp(2) min(df)];
    case 'multiple linear'
        len = length(data(:,1));
        a = cell(ops,5);
        for i = 1:ops
            a{i,5} = regressors;
            a{i,1} = regress(data(:,i),a{i,5});
            a{i,2} = 1;
            a{i,3} = 1;
            a{i,4} = (a{i,5}'*a{i,5})./(len.^2);
        end
    case 'multivariate normal'
        MV = robustmcov(data);
        a = cell(1,4);
        a{1} = MV{1};
        a{2} = 1;
        a{3} = 1;
        a{4} = MV{2};
end
end
