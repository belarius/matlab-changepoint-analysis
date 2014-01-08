function P = EstimateParameters(data,distr,M,knownp,a,PARAM)
%ESTIMATEPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
ops = size(data,2);
switch distr
    case {'linear','multiple linear','multivariate normal'}
        P = cell(length(M)-1,2);
    case {'binomial','geometric'}
        P = zeros(length(M)-1,1);
    otherwise
        P = zeros(length(M)-1,ops);
end
for i = 1:length(M)-1
    d = data(M(i)+1:M(i+1),:);
    len = M(i+1)-M(i);
    ops = length(data(1,:));
    smd = sum(d);
    if PARAM==0 %==BAYESIAN PARAMETER ESTIMATES==
        q = PosteriorHyper(d,distr,a,knownp,[M(i)+1 M(i+1)]);
        switch distr
            case 'binomial'
                P(i,1) = (a(1) + smd(1))./(sum(a) + sum(smd));
            case 'geometric'
                P(i,1) = (a(1) + len)./(a(1) + a(2) + len + smd);
            case 'poisson'
                P(i,1) = (a(1) + smd)./(a(2) + len);
            case 'multinomial'
                P(i,:) = (sum(d,1)+a)./(sum(a) + sum(smd));
            case 'exponential'
                P(i,1) = (a(2) + smd)./(a(1) + len);
            case 'linear'
                P{i,1} = q{1};
                P{i,2} = inv(q{4});
                f = a{5}(M(i)+1:M(i+1),:)*P{i,1};
                P{i,2} = (std(d-f).^2).*P{i,2};
            case 'normal'
                if isnan(knownp(1)) && isnan(knownp(2))
                    P(i,1) = q(1);
                    P(i,2) = sqrt((q(4).*(q(2)+1))./(q(2).*q(3)));
                elseif isnan(knownp(1))
                    prec = 1./(knownp(2).^2);
                    P(i,1) = (a(1).*a(2) + len.*prec.*mean(d))./(a(2) + len.*prec);
                    P(i,2) = knownp(2);
                elseif isnan(knownp(2))
                    P(i,1) = knownp(1);
                    P(i,2) = sqrt((a(2) + sum((d - knownp(1)).^2)./2)./(a(1) + len./2));
                end
            case 'uniform'
                if isnan(knownp(1))
                    P(i,1) = min([a(2);d]);
                else
                    P(i,1) = knownp(1);
                end
                if isnan(knownp(2))
                    P(i,2) = max([a(3);d]);
                else
                    P(i,2) = knownp(2);
                end
            case 'multiple linear'
                P{i,1} = cell(1,ops);
                P{i,2} = cell(1,ops);
                for j = 1:ops
                    if M(i)+1 == M(i+1)
                        P{i,1}(j) = {q{j,1}};                        
                        P{i,2}(j) = {nan};
                    else
                        P{i,1}(j) = {q{j,1}};
                        warning('off','all');
                        P{i,2}(j) = {inv(q{j,4})};
                        warning('on','all');
                        f = a{j,5}(M(i)+1:M(i+1),:)*P{i,1}{j};
                        P{i,2}{j} = (std(d(:,j)-f).^2).*P{i,2}{j};
                        trash = 0;
                    end
                end
            case 'multivariate normal'
                P{i,1} = q{1};
                P{i,2} = q{4}./(q{3} - ops - 1);
        end
    elseif PARAM==1 %==FREQUENTIST PARAMETER ESTIMATES==
        switch distr
            case {'binomial','poisson','exponential'}
                P(i,1) = mean(d(:,1));
            case {'multinomial'}
                P(i,:) = mean(d);
            case 'geometric'
                P(i,1) = len./(len + smd);
            case 'linear'
                [P{i,1},~,~,P{i,2}] = mvregress(a{5}(M(i)+1:M(i+1),:),d);
            case 'normal'
                if isnan(knownp(1)) && isnan(knownp(2))
                    P(i,1) = mean(d);
                    P(i,2) = std(d);
                elseif isnan(knownp(1))
                    P(i,1) = knownp(1);
                    P(i,2) = std(d);
                elseif isnan(knownp(2))
                    P(i,1) = mean(d);
                    P(i,2) = knownp(2);
                end
            case 'uniform'
                if isnan(knownp(1))
                    P(i,1) = min(d);
                else
                    P(i,1) = knownp(1);
                end
                if isnan(knownp(2))
                    P(i,2) = max(d);
                else
                    P(i,2) = knownp(2);
                end
            case 'multiple linear'
                for j = 1:ops
                    [P{i,1}(j),~,~,P{i,2}(j)] = mvregress(a{5}(M(i)+1:M(i+1),:),d(:,j));
                end
            case 'multivariate normal'
                P{i,1} = mean(d);
                P{i,2} = cov(d);
        end
    elseif PARAM==2 %==ROBUST PARAMETER ESTIMATES==
        switch distr
            case {'binomial'}
                P(i,1) = mean(d(:,1));
            case {'multinomial'}
                P(i,:) = mean(d);
            case 'geometric'
                P(i,1) = len./(len + smd);
            case {'poisson','exponential'}
                P(i,1) = median(d);
            case 'linear'
                [P{i,1},P{i,2}] = robustfit(a{5}(M(i)+1:M(i+1),:),d,'bisquare',4.685,'off');
                P{i,2} = P{i,2}.covb;
            case 'normal'
                if isnan(knownp(1)) && isnan(knownp(2))
                    P(i,1) = median(d);
                    P(i,2) = mad(d,1)*1.4826;
                elseif isnan(knownp(1))
                    P(i,1) = knownp(1);
                    P(i,2) = mad(d,1)*1.4826;
                elseif isnan(knownp(2))
                    P(i,1) = median(d);
                    P(i,2) = knownp(2);
                end
            case 'uniform'
                if isnan(knownp(1))
                    P(i,1) = min(d);
                else
                    P(i,1) = knownp(1);
                end
                if isnan(knownp(2))
                    P(i,2) = max(d);
                else
                    P(i,2) = knownp(2);
                end
            case 'multiple linear'
                for j = 1:ops
                    [P{i,1}(j),~,~,P{i,2}(j)] = mvregress(a{5}(M(i)+1:M(i+1),:),d(:,j));
                    [P{i,1}(j),P{i,2}(j)] = robustfit(a{5}(M(i)+1:M(i+1),:),d(:,j),'bisquare',4.685,'off');
                    P{i,2}(j) = P{i,2}(j).covb;
                end
            case 'multivariate normal'
                MV = robustmcov(d);
                P{i,1} = MV{1};
                P{i,2} = MV{2};
        end
    end
end
end
