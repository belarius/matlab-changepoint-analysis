function q = PosteriorHyper(data,distr,a,knownp,range)
%POSTERIORHYPER Summary of this function goes here
%   Detailed explanation goes here
switch distr
    case 'uniform'
        if ~isnan(knownp(1))
            data = data-knownp(1);
            a(2) = 0;
            a(3) = a(3)-knownp(1);
            a(4) = a(3);
        elseif ~isnan(knownp(2))
            data = knownp(2)-data;
            a(3) = knownp(2)-a(2);
            a(2) = 0;
            a(4) = a(3);
        elseif ~isnan(a(3)-a(2))
            a(4) = a(3)-a(2);
        end
end

%posterior hyperparameters
len = length(data(:,1));
switch distr
    case 'binomial'
        q = sum(data)+a;
    case 'geometric'
        q = [a(1)+len a(2)+sum(data)];
    case 'poisson'
        q = [a(1)+sum(data) a(2)+len];
    case 'multinomial'
        q = sum(data)+a;
    case 'exponential'
        q = [a(1)+len a(2)+sum(data)];
    case 'linear'
        q = local_linear_posth(data,len,a,range);
    case 'normal'
        if isnan(knownp(1)) && isnan(knownp(2))
            q = local_normal_mv_posth(data,len,a);
        elseif isnan(knownp(1))
            q = local_normal_m_posth(data,len,a,knownp);
        elseif isnan(knownp(2))
            q = local_normal_v_posth(data,len,a,knownp);
        end
    case 'uniform'
        q = local_uni_posth(data,len,a,knownp);
    case 'multiple linear'
        q = local_multilinear_posth(data,len,a,range);
    case 'multivariate normal'
        q = local_normal_invwishart_posth(data,len,a);
end
    
end

%====POSTERIOR HYPERPARAMETER SUBROUTINES====

function q = local_normal_m_posth(d,len,a,knownp)
    q = zeros(1,2);
    q(2) = (a(2) + len)./(knownp(2).^2);
    q(1) = (a(1)./(knownp(2).^2) + sum(d)./(knownp(2).^2)) ./ q(2);
end

function q = local_normal_v_posth(d,len,a,knownp)
    q = zeros(1,2);
    q(1) = a(1) + len./2;
    q(2) = a(2) + sum((d-knownp(1)).^2)./2;
end

function q = local_normal_mv_posth(d,len,a)
    q = zeros(1,4);
    mn = mean(d);
    q(1) = (a(2).*a(1) + len.*mn)./(a(2) + len);
    q(2) = a(2) + len;
    q(3) = a(3) + len./2;
    q(4) = a(4) + 0.5.*sum((d-mn).^2) + ((len.*a(2).*((mn - a(1)).^2))./(2.*(a(2) + len)));
end

function q = local_linear_posth(d,len,a,range)
    warning('off','MATLAB:singularMatrix')
    x = a{5}(range(1):range(2),:);
    q = cell(1,4);
	q{4} = x'*x + a{4};
	q{2} = a{2} + len./2;
    q{1} = (x'*x + a{4})\(x'*d + a{4}*a{1}); % (x'*x + a{4})\(x'*d + a{4}*a{1}) is equivalent to inv(x'*x + a{4})*(x'*d + a{4}*a{1})
	q{3} = a{3} + 0.5.*(d'*d + a{1}'*a{4}*a{1} - q{1}'*q{4}*q{1});
    warning('on','MATLAB:singularMatrix')
end

function q = local_multilinear_posth(d,len,a,range)
    ops = length(d(1,:));
    q = cell(ops,4);
    for i = 1:ops
        q(i,:) = local_linear_posth(d(:,i),len,a(i,:),range);
    end
end

function q = local_uni_posth(d,len,a,knownp)
    if length(d)==1
        q =  nan(1,4);
    else
        q(1) = a(1) + len;
        if ~isnan(a(2)) && ~isnan(a(3))
            q(2) = min([d;a(2)]);
            q(3) = max([d;a(3)]);
            q(4) = q(3)-q(2);
        else
            q(2) = NaN;
            q(3) = NaN;
            q(4) = max([a(4);max(d)-min(d)]);
        end
    end
end

function q = local_normal_invwishart_posth(d,len,a)
    q = cell(1,4);
    q{1} = (a{2}./(a{2}+len))*a{1} + (len./(a{2}+len))*mean(d);
    q{2} = a{2} + len;
    q{3} = a{3} + len;
    q{4} = inv(a{4}) + cov(d).*len + ((len.*a{2})./(len + a{2}))*(mean(d)-a{1})*(mean(d)-a{1})';

end

