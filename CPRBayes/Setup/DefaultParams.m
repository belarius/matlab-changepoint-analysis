function knownp = DefaultParams(distr,ops)
%DEFAULTPARAMS Sets the known parameters for the model in the event they aren't set by the user.
%   ~fill in later~
switch distr
    case {'binomial','geometric','poisson','exponential'}
        knownp = NaN;
    case {'linear','normal','uniform'}
        knownp = [NaN NaN];
    case {'multiple linear','multivariate normal'}
        knownp = cell(1,2);
        knownp{1} = NaN;
        knownp{2} = NaN;
    case 'multinomial'
        knownp = nan(1,ops);
    otherwise
        error('CPRBayes:UnknownDistr','Distribution not recognized');
end
end
