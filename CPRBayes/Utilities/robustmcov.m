function MV = robustmcov(x)
%ROBUSTMCOV Robust estimation of the mean and covariance, (Campbell, 1980)
%   Campbell, N. A. (1980) Robust procedures in multivariate analysis I:
%       Robust covariance estimation. Journal of the Royal Statistical
%       Society, Series C (Applied Statistics), 29(3), 231-237.

MV = cell(1,2);
[n v] = size(x); xm = mean(x);
xdiff = x - repmat(xm,n,1);
b1 = 2; b2 = 1.25;
d0 = sqrt(v) + b1./sqrt(2);
w = zeros(n,1);
for i = 1:n
    dm = sqrt( (xdiff(i,:)')'*inv(cov(x))*(xdiff(i,:)') );
    if dm <= d0
        w(i) = 1;
    else
        w(i) = (d0.*exp(-(dm-d0).^2/(2.*b2.^2)))./dm;
    end
end
num1 = zeros(1,v); denom1=0;
for i = 1:n
    num1 = num1 + w(i).*x(i,:);
    denom1 = denom1 + w(i);
end
num2 = zeros(v); denom2=0;
MV{1} = num1./denom1;
xdiff = x - repmat(MV{1},n,1);
for i = 1:n
    num2 = num2 + w(i).^2 * (xdiff(i,:)')*(xdiff(i,:)')';
    denom2 = denom2 + w(i).^2;
end
MV{2} = num2./(denom2-1);
end

