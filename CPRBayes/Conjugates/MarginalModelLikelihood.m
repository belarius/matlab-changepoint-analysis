function MML = MarginalModelLikelihood(data,distr,q,a,knownp,ops)
%MARGINALMODELLIKELIHOOD Calculate MML based on the appropriate conjugate prior
%   Detailed explanation goes here
switch distr
    case {'binomial','geometric','multinomial'}
        MML = local_multinom_nc_ln(q);
    case {'poisson','exponential'}
        MML = local_gamma_nc_ln(q);
    case 'linear'
        if length(data) < 3
            MML = -inf;
        else
            MML = local_linear_nc_ln(q,a);
        end
    case 'normal'
        if isnan(knownp(1)) && isnan(knownp(2))
            MML = local_normalgamma_nc_ln(q);
        elseif isnan(knownp(1))
            MML = local_normal_nc_ln(q);
        elseif isnan(knownp(2))
            MML = local_gamma_nc_ln(q);
        end
    case 'uniform'
        if isnan(knownp(1)) && isnan(knownp(2))
            MML = local_uni_d_ln(q);
        else
            MML = local_uni_s_ln(q);
        end
    case 'multiple linear'
        MML = local_multilinear_nc_ln(q,a);
    case 'multivariate normal'
        if size(data,1) >= ops
            MML = local_normal_invwishart_nc_ln(ops,q);
        elseif size(data,1) >= 1
            m = q{1};
            v = a{4};%0.5*eye(ops) + 0.5*ones(ops);
            MML = -(ops./2).*log(ops);
            for i = 1:size(data,1)
                MML = MML - (ops./2).*log(2.*pi) - 0.5.*log(det(v));
                MML = MML - 0.5.*(data(i,:)-m)*inv(v)*(data(i,:)-m)';
            end
        else
            v = a{4};%0.5*eye(ops) + 0.5*ones(ops);
            MML = (-ops./2).*log(2.*pi) - 0.5.*log(det(v));
            MML = MML + (-ops./2).*log(ops);
        end
end
end

%====NORMALIZING CONSTANT SUNROUTINES====
function out = local_gamma_nc_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of the gamma
%   a(1) = shape
%   a(2) = rate
    out = gammaln(a(1)) - a(1).*log(a(2));
end

function out = local_uni_s_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of a uniform
%distribution with one anchored tail.
%   a  = parameters
    out = -1.*( log(a(1)) + a(1).*log(a(2)));
end

function out = local_uni_d_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of the uniform
%distribution with neither tail anchored
%   a  = parameters
    out = -1.*( log(a(1)) + log(a(1)+1) + a(1).*log(a(4)));
end

function out = local_multinom_nc_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of the dirichlet-multinomial
%   a  = parameters
    if ~isempty(find(a==0, 1))
        out = 0;
    else
        out = sum(gammaln(a)) - gammaln(sum(a));
    end
end

function out = local_normal_nc_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of the normal
%   a  = parameters
    out = log(sqrt(2.*pi./a(2))) + (a(1).^2).*a(2)./2;
    %vare = (stdp.^2)./a(2);
    %out = 0.5.*log(2.*pi.*vare) + (a(1).^2)./(2.*vare);
end

function out = local_normalgamma_nc_ln(a)
%local_multinom_nc_ln Returns the log normalizing constant of the normal-gamma
%   a  = parameters
    out = gammaln(a(3)) - a(3).*log(a(4)) + 0.5.*(log(2.*pi) - log(a(2)));
end

function out = local_linear_nc_ln(q,a)
    len = 2.*(q{2} - a{2});
    out = (-len./2).*log(2.*pi()) + 0.5.*log(det(a{4})./det(q{4})) + a{2}.*log(a{3}) - q{2}.*log(q{3}) + gammaln(q{2}) - gammaln(a{2});
end

function out = local_multilinear_nc_ln(q,a)
    out = 0;
    for i = 1:length(q(:,1))
        qx = q(i,:);
        ax = a(i,:);
        out = out + local_linear_nc_ln(qx,ax);
    end
end

function out = local_normal_invwishart_nc_ln(ops,a)
%local_multinom_nc_ln Returns the log normalizing constant of the gamma
%   a  = parameters

out = 0;
out = out + 0.5.*a{3}.*ops.*log(2);
out = out + 0.25.*ops.*(ops-1).*log(pi());
for i = 1:ops
    out = out + gammaln((a{3} + 1 - i)./2);
end
out = out + 0.5.*ops.*log((2.*pi())./a{2});
out = out - 0.5.*a{2}.*log(det(a{4}));

end

