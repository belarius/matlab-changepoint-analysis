function [V] = CPVerbose(M,R)
%CPVerbose converts the compact form of a binomial changepoint analysis to a verbose form.
%   Details comment

len = M(end);
V = zeros(1,len);
dex = 2;
for i = 1:len
    if M(dex) < i
        dex = dex+1;
    end
    V(i) = R(dex-1);
end

