function [M,P,stats] = CPRBayes(D,distr,varargin)
%CPRBAYES Estimates a changepoint model using marignal likelihood ratios.
%   [M,P] = CPRBayes(D,DISTR) performs model comparisons to divide the time
%   series D into discontinuous subsections, divided by changepoints reported
%   in M, fit to distribution DISTR, and represented by parameters P.
%
%   D is a matrix with rows corresponding to observations, and DISTR is the
%   distribution believed to be responsible for generating observations D.
%   The function interprets D differently depending on this distribution. D
%   is assumed to represent a uniform subdivision of time/trials.
%   Observations in D are also assumed to be exchangeable within each
%   changepoint-delimited subdivision.
%
%   DISTR must be specified by the user, and may be any of the following:
%     ==Discrete Distributions==
%       'binomial' - The data in D are presumed to *either* be individual
%               observations of a Bernoulli process, consisting of the
%               values 0 and 1 in a single column, *or* to consist of 
%               binomial data with k successes in data(:,1) out of n trials
%               in data(:,2).
%       'geometric' - The data in D are presumed to represent the number of
%               failures before the first success on a series of Bernoulli
%               trials. D should have a single column of integer data.
%       'poisson' - The data in D are presumed to be generated by a Poisson
%               process. D should have a single column of integer data.
%       'multinomial' - The multivariate generalization of the binomial
%               distribution. If D has a single column, it is presumed to
%               be categorical, with each unique value in D constituting a
%               number of observation in a respective category. If D has
%               more than one column, each column is presumed to be a
%               logical vector corresponding to a predefined category.
%     ==Continuous Distributions==
%       'exponential' - The data in D are presumed to be drawn from an
%               exponential distribution. D should have a single column.
%       'linear' - The data D are presumed to result in an additive fashion
%               from a set of linear predictors. D should have a single
%               column.
%       'normal' - The data in D are presumed to be drawn from a Gaussian
%               distribution whose means and precision may either be
%               specified, but are unknown by default. D should have a
%               single column.
%       'uniform' - The data in D are presumed to be drawn from a uniform
%               distribution.
%       'multivariate normal' - The data in D are presumed to be drawn from
%               a multivariate Gaussian distribution
%       'multiple linear' - The data in D are presumed to be a series of
%               independent outcomes all predicted by a common set of
%               explanatory linear factors. This is *not* a General Linear
%               Model, as it lacks regression parameters for full
%               covariance. D should consists of 2 or more columns, each of
%               which is an outcome being modeled independently.
%
%   [M,P,stats] = CPRBayes(D,DISTR,'PARAM1',val1,'PARAM2',val2,...) allows you to
%   specify optional parameter name/value pairs to fix certain aspects of
%   the analysis, rather than letting them be inferred from default
%   settings. Parameters are:
%
%      'alpha' - specify the array of hyperparameters that will be used to
%           estimate the normalizing constant in the ratio of marginal
%           model likelihoods. Accepts a row vector of parameter values.
%           See "parameter specifications" below for details associated
%           with a particular distribution.
%
%      'conf' - specify the confidence interval. Default is .95,
%           corresponding to the 95% interval.
%
%      'depth' - Default value of zero, which results in normal behavior.
%           Any higher value forces the algorithm to perform a fixed number
%           of subdivisions, regardless of whether the evidence supports
%           additional cps. This can be useful when trying to identify the
%           `least worst' cps, even though the evidence does not favor
%           them.
%
%      'knownp' - specify parameters whose value is known. While not
%           required for the supported distributions, having known
%           parameters changes the marginal solution. At present, this
%           setting only applies to the normal and uniform models.
%
%      'knownM' - specify known subdivisions in the data. The specified
%           model M is expected to be a vector with at least values
%           values [0 len], as well as intermediate specified points.
%
%      'priorc' - specify the prior probability of a change per
%           observation. If left unspecified, the prior probability is
%           assumed to be 1./(length(D)-1) at the outset, and grows as
%           changepoints are identified. When specified with this
%           parameter, the probability is instead held constant at the
%           given value.
%
%      'timestamps' - specify an array of timestamps corresponding to the
%           observations. Times are presumed to appear in a uniform fashion
%           by default.
%
%      'regressors' - applies only to the linear and multiple linear models
%           (the input is replaced by NaN otherwise). This specifies the
%           regression's explanatoryt factors. By default, a model with an
%           intercept and the timestamps variable is used.
%
%      'thresh' - specify the threshold criterion used to assess whether a
%           changepoint is likely given the evidence. By default, any value
%           greater than 10 is taken as evidence in favor of a changepoint.
%           Since the threshold is compared to the posterior odds of a
%           changepoint being present, it should be >= 1; furthermore, the
%           algorithm performs poorly if a threshold less than 3 is used.
%
%      'weight' - specify a weight for each candidate interval, which takes
%           the form ones(len,1). These values are multiplied by the
%           marginal likelihoods. If an interval is known to be impossible,
%           its weight may be set to zero. Weights greater than 1 should be
%           avoided, as they risk compromising the conservatism of the
%           analysis.
%
%      'messages' - If 'on' (default), then post messages; otherwise, no
%           messages.
%
%      'conftype' - Accepts one of two strings: 'itemwise' (default)
%           or 'incremental'. 'incremental' uses the CIs that emerge as a
%           result of the pairwise process; 'itemwise' is a post-hoc
%           approach that assesses the CI on the basis of the bookending
%           CPs once all CPs have been identified. As a rule of thumb, the
%           'itemwise' CIs will overlap less, since they take the existence
%           of the CPs as a given; by contrast, the 'incremental' CIs use
%           larger ranges of data, but are likely to show dramatic and
%           irregular overlap. 'itemwise' is also in principle more
%           consistent with the computed marginal likelihoods
%
%      'paramtype' - Accepts one of three strings: 'bayesian' (default),
%           'frequentist', or 'robust'. 'bayesian' follows the logic of
%           Bayesian updating to make posterior parameter estimates that
%           incorporate the prior. 'frequentist' makes use of traditional
%           parameter estimation techniques, such as arithmetic averaging.
%           'robust' makes use of robust parameter estimation techniques.
%
%      'dice' - An exploratory technique used if the data are suspected to
%           consist of smaller sub-processes embedded within larger
%           stationary or cyclical processes. If dice is invoked, the
%           integer value it is given indicates the number of equally-sized
%           segments into which the data are divided. These segments are
%           then checked for a single change-point each, and any
%           change-points that are discovered are retained. The algorithm
%           then proceeds as normal, with knownM reflecting the additional
%           change-points.
%
%===OVERALL LOGIC==
%
%   See Jensen (Submitted) for details.
%
%===MODEL INTERPRETATION==
%
%   M is a column vector listing the last items in each subsection. In all
%   cases, M(1) = 0, and each subdivision of the data spans M(i)+1:M(i+1).
%   Changepoints are presumed to reside between M(i) and M(i)+1, and the
%   vector denotes M(i) as the changepoint merely as a naming convention.
%
%   P is a matrix where each row corresponds to the model parameters in a
%   particular subdivision of the data, and each column refers to a
%   specific distribution parameter for that segment. Given a matrix
%   P(i,j), the identity of each column j varies, as specified below:
%
%     ==Discrete Distribution P parameters==
%       'binomial'
%           P(i,1) = p (Probability of success)
%       'geometric'
%           P(i,1) = p (Probability of success)
%       'poisson'
%           P(i,1) = lambda (rate)
%       'multinomial'
%           P(i,j) = p (Probability of category j), where sum(P(i,:)) = 1
%     ==Continuous Distribution as==
%       'exponential'
%           P(i,1) = lambda (mean)
%       'normal'
%           P(i,1) = mu (mean)
%           P(i,2) = sigma (standard deviation)
%       'linear'
%           P{i,1} = beta (vector of regression coefficients)
%           P{i,2} = covariance error matrix
%       'uniform'
%           P(i,1) = mn (minimum)
%           P(i,2) = mx (maximum)
%       'multiple linear'
%           P{i,1}{j} = beta for outcome j (vector of regression coefficients)
%           P{i,2}{j} = covariance error matrix for outcome j
%       'multivariate normal'
%           P{i,1} = Mu (vector of mean)
%           P{i,2} = Sigma (covariance matrix)
%
%===HYPERPARAMETER SPECIFICATION==
%
%   In order to infer the noramlizing constant of a distribition's
%   conjugate prior, it is necessary to estimate the appropriate posteiror
%   hyperparameters. These are defined in part by the evidence, but also by
%   the analyst's prior assumptions. In this function, we denote the prior
%   hyperparameters using the vector a (short for 'alpha'), which contains
%   two or more parameters, whose definition differs for each distribution.
%
%   These parameters can be considered "virtual observations" in some
%   cases, to be intermixed with actual observations. By default, all
%   distributions begin with all values that constitutes weak priors.
%   For example, in the binomial case, a(1) = a(2) = 0.5 is a sufficiently
%   weak prior that it may be considered a reference prior. However, if an
%   analyst has good reasons to think the 50/50 odds are correct, then a
%   prior of a(1) = a(2) = 100 would constitute a very strong expectation
%   of equality, effectively adding 200 observations to the data.
%   
%   In most cases, an analyst should rely on either strong theoretical or
%   empirical reasons when using strong prior hyperparameters, as these
%   increase the odds of washing out small but real transitions in
%   behavior.
%
%   In the list below, we indicate both the interpretation of the prior
%   hyperparameters and the equation used to estimate the posterior
%   hyperparameters, which are denoted as q(). Other shorthand
%   includes n = length(D); xbar = mean(D);
%   
%   Unless otherwise noted, support for these is provided by DeGroot
%   (2004).
%   
%     ==Discrete Distribution as==
%       'binomial' = [0.5 0.5]
%           Prior Hyperparamters
%               a(1) = ("a" - 1) successes
%               a(2) = ("beta" - 1) failures
%           Posterior Hyperparameters
%               q(1) = a(1) + sum(D)
%               q(2) = a(2) + n - sum(D)
%           Normalizing Constant
%               NC = gamma(a(1) + a(2))./(gamma(a(1)).*gamma(a(2))
%           Posterior Predictive
%               p = q(1)./(q(1)+q(2))
%       'geometric' = 0.5 0.5
%           Prior Hyperparamters
%               a(1) = ("a" - 1) experiments
%               a(2) = ("beta" - 1) failures
%           Posterior Hyperparameters
%               q(1) = a(1) + n
%               q(2) = a(2) + sum(D)
%           Normalizing Constant
%               NC = gamma(a(1) + a(2))./(gamma(a(1)).*gamma(a(2))
%           Posterior Predictive
%               p = (q(1)+n)./(q(1)+q(2)+sum(D)+n)
%       'poisson' = [mean(D) 1]
%           Prior Hyperparamters
%               a(1) = "a" total occurances
%               a(2) = "beta" intervals
%           Posterior Hyperparameters
%               q(1) = a(1) + sum(D)
%               q(2) = 1./(a(2) + n)
%           Normalizing Constant
%               NC = (beta.^a).*gamma(a)
%           Posterior Predictive
%               lambda = ???
%       'multinomial' = 0.5.*ones(i,1)
%           Prior Hyperparamters
%               a(i) = ("a_i" - 1) occurances of category i
%           Posterior Hyperparameters
%               q(i) = a(i) + sum(D(:,i))
%           Normalizing Constant
%               NC = gamma(sum(a))./product(gamma(a))
%           Posterior Predictive
%               P = dirichlet(q) 
%     ==Continuous Distribution as==
%       'exponential' = [1 median(D)]
%           Prior Hyperparamters
%               a(1) = "a" observations
%               a(2) = "beta" sum of observations
%           Posterior Hyperparameters
%               q(1) = a(1) + n
%               q(2) = a(2) + sum(D)
%           Normalizing Constant
%               NC = (beta.^a).*gamma(a)
%           Posterior Predictive
%               lambda = (q(2))./(q(1))
%       'normal' = [median(D) 0.1 0.1 mad(D,1).*1.4826]
%           Prior Hyperparamters
%               a(1) = "mu_0" population mean; is estimated robustly using the median
%               a(2) = "v" observations for the mean
%               a(3) = (2.*"a" + 1) observations for the precision
%               a(4) = "beta" sum of squares; is estimated robustly using the median absolute deviation
%           Posterior Hyperparameters
%               q(1) = (a(1).*a(2) + n.*xbar)./(a(2) + n)
%               q(2) = a(2) + n
%               q(3) = a(3) + (n./2)
%               q(4) = a(4) + 0.5.*sum((D-xbar).^2) +
%                           (n.*a(2))./(a(2) + n) .*
%                           ((xbar - a(1)).^2)./2
%           Normalizing Constant
%               NC = (gamma(a)./(beta.^a)).*sqrt((2.*pi)./v)
%           Posterior Predictive
%               mu = q(1);
%               sigma = sqrt((q(4).*(q(2)+1))./(q(2).*q(3)))
%       'uniform' = [1 knownMin knownMax knownGap]
%           Prior Hyperparamters
%               a(1) = "a" observations
%               a(2) = "r1" minimum; not jointly specified with a(3)
%               a(3) = "r2" maximum; not jointly specified with a(2)
%               a(4) = "r2-r1" range; can be specified independently of a(2) and a(3); is arbitrarily set to 1.0 if no parameters are known
%           Posterior Hyperparameters
%               q(1) = a(1) + n
%               q(2) = min([a(2);D])
%               q(3) = max([a(3);D])
%           Normalizing Constant
%               NC = a(1).*(a(1)+1).*((a(3)-a(2)).^a(1))
%           Posterior Predictive
%               mn = min([a(2);D])
%               mx = max([a(3);D])
%       'multivariate normal' = {mean_hat 1 1 covariance_hat}; these are estimated using a robust estimation procedure described by Campbell (1980). 
%           Prior Hyperparamters
%               a{1} = mean vector
%               a{2} = "k" observations supporting the mean 
%               a{3} = "v" observations supporting the covariance
%               a{4} = covariance matrix
%           Posterior Hyperparameters
%               q{1} = (a{2}./(a{2}+n))*a{1} + (n./(a{2}+n))*xbar;
%               q{2} = a{2} + n;
%               q{3} = a{3} + n;
%               q{4} = inv(a{4}) + cov(D).*n + ((n.*a{2})./(n + a{2}))*(xbar-a{1})*(xbar-a{1})';
%           Normalizing Constant & Posterior Predictive
%               These have a lot of terms. See the treatment of the
%               Normal-Inverse-Wishart in DeGroot (2004).
%       'linear' and 'multiple linear'
%           The implementation of Bayesian linear regression vis conjugate
%           priors is rather involved. See Gelman et al. (2003) for
%           details.
%
%===REFERENCES===
%   Campbell NA (1980) Robust procedures in multivariate analysis I:
%       Robust covariance estimation. Journal of the Royal Statistical
%       Society, Series C (Applied Statistics), 29(3), 231-237.
%   DeGroot MH (2004) Optimal Statistical Decisions. Wiley Classics
%       Library.
%   Gelman A, Carlin JB, Stern HS, Rubin DB (2003) Bayesian Data Analysis,
%       2nd edition. CRC PRess.
%   Jensen G (Submitted) Closed-form estimation of multiple change-points
%       models.
%
% written by:
% Greg Jensen
% Columbia University (Psychology)
% greg.guichard.jensen@gmail.com

%===INITIALIZATION===
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

data = SetupData(D,distr);
SETTINGS = SetupSettings(data,distr,varargin);

nanv = find(sum(isnan(data),2)==1);
if ~isempty(nanv)
    data(nanv,:) = [];
end

len = SETTINGS.LEN;
ops = SETTINGS.OPS;
p_c = SETTINGS.PRIOR_C;
M = SETTINGS.PRIOR_M;
times = SETTINGS.TIME;
regressors = SETTINGS.REGRESSORS;
thresh = SETTINGS.THRESH;
weight = SETTINGS.WEIGHT;

if SETTINGS.CONF
    Conf_Placeholder = SETTINGS.CONF;
    SETTINGS.CONF = 0;
end

ticID = tic;

R = zeros(length(M),1);

%===FINAL SETUP ERROR CHECKING===

%=is there even data?=
if len == 0
    error('CPRBayes:DataLength','Data variable D has 0 rows');
end

%=size errors=
if size(times,1) ~= len
    error('CPRBayes:TimestampDataLengthMismatch','times() has %i timestamps when it should have %i timestamps',length(times),len);
end

if strcmp(distr,'linear')||strcmp(distr,'multiple linear')
    if size(regressors,1) ~= len
        error('CPRBayes:TimestampDataLengthMismatch','regressors() has %i values when it should have %i values',length(regressors),len);
    end
end
if size(times,2) > 1
	error('CPRBayes:TimestampsMustBeUnivariate','times() has %i colums but should only have 1 column',size(times,2));
end

%===DICING===

if SETTINGS.DICE > 1
    [M,R] = CPRBayes_dice(data,distr,SETTINGS);
end

%===MAIN RECURSIVE ALGORITHM==
count = 1;                      %Denotes the number of additional cps found in each cycle. Begins at 1 to initiate the loop.
deep = 0;
k = [];
td= [];
q = ones(len,1);
while count > 0;                %When no new cps are identified, the loop ends.
    m = nan(length(M)-1,1);   %Denotes the different subsections of the data, each of which is checked for additional cps.
    r = nan(length(M)-1,1);
    deep = deep+1;
    k(:,deep) = nan(len-1,1);
    td(:,deep) = nan(len-1,1);
    q(:,deep+1) = ones(len,1);
    for i = 1:length(M)-1
        if q(M(i)+1,deep) == 1
            segment = data(M(i)+1:M(i+1),:);
            segtime = times(M(i)+1:M(i+1),:);
            segweight = weight(M(i)+1:M(i+1)-1,:);
            [cp,~,k_c,td_c] = CPRBayes_single(segment,distr,SETTINGS,segtime,segweight);   %The identity and Bayes Factor of a putative cp is identified in segment i.
            BF = exp(k_c + td_c).*segweight;
            BFs = sum(BF(~isnan(BF)));
            k(M(i)+1:M(i+1)-1,deep) = k_c;
            td(M(i)+1:M(i+1)-1,deep) = td_c;
            if ~isnan(BFs) && ~isnan(cp)
                p = p_c(1)./p_c(2);
                posterior = BFs.*p.*((M(i+1))-(M(i)+1));
                m(i) = cp;
                r(i) = posterior;
            end
        else
            k_c = k(M(i)+1:M(i+1)-1,deep-1);
            k(M(i)+1:M(i+1)-1,deep) = k_c;
            td_c = td(M(i)+1:M(i+1)-1,deep-1);
            td(M(i)+1:M(i+1)-1,deep) = td_c;
            BF = exp(k_c + td_c).*weight(M(i)+1:M(i+1)-1,:);
            BFs = sum(BF(~isnan(BF)));
            k(M(i)+1:M(i+1)-1,deep) = k_c;
            td(M(i)+1:M(i+1)-1,deep) = td_c;
            if ~isnan(BFs) && ~isempty(BF)
                p = p_c(1)./p_c(2);
                posterior = BFs.*p.*(times(M(i+1))-times(M(i)+1));
                m(i) = find(BF==max(BF),1);
                r(i) = posterior;
            end
        end
    end
    mtest = ~isnan(m);        % Where were the MLE changes?
    rtest = r>thresh;   % Where were the MLE changes whose MML justifies their inclusion?
    if deep > 1
        test2 = sum(k(:,deep)~=k(:,deep-1));
    else
        test2 = 1;
    end
    if ~test2
        count = 0;
    else
        count = length(find(rtest));
    end
    qm = m + M(1:length(M)-1);      %This step updates the cp indices to reflect the global data indexing, rather than local positions in each segment.
    if ~SETTINGS.DEPTH
        for i = find(~rtest)
            q(M(i)+1:M(i+1),deep+1) = 0;
        end
        if count>0
            sortblock = sortrows([M R;qm(rtest) r(rtest)],1);
        else
            sortblock = sortrows([M R],1);
        end
    else
        if sum(mtest)>0
            sortblock = sortrows([M R;qm(mtest) r(mtest)],1);
        else
            sortblock = sortrows([M R],1);
        end
    end
    M = sortblock(:,1);
    R = sortblock(:,2);
    mlist = unique(M);
    if sum(histc(M,mlist)>1)>0
        for i = 1:length(mlist)
            del = find(M==mlist(i));
            del(R(del)==max(R(del))) = [];
            M(del) = [];
            R(del) = [];
        end
    end

    p_c(1) = length(M)-2;       %When p_c is not fixed, the prior odds of finding a cp is updated at the end of each loop.

    if SETTINGS.DEPTH
        if deep == SETTINGS.DEPTH
            count = 0;
        else
            count = 1;
        end
    else
        if count == 0
            deep=deep-1;
        end
    end
end

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

stats = struct('model',{M},'distr',{distr},'params',{P},'timestamps',times,'regressors',regressors,'weights',weight,'threshold',{thresh},'post_ratios',{R},'seg_numer',k,'seg_denom',td,'conf',SETTINGS.CONF_A,'conf_int',{C},'prior_change',{p_c(1)./p_c(2)},'prior_model',{SETTINGS.ALPHA},'depth',{deep},'elapsed_time',ElapsedTime);

end