function [cp,ci,k_c,td_c] = CPRBayes_single(data,distr,SETTINGS,times,weight)
%CPRBAYES_SINGLE Summary of this function goes here
%   Detailed explanation goes here

a = SETTINGS.ALPHA;
ops = SETTINGS.OPS;
knownp = SETTINGS.KNOWN_P;
conf = SETTINGS.CONF_A;

len = length(data(:,1));
%==Calculate bayes factors==
if len == 1 %since CPs can only happen between points, the odds given no such interval are undefined
    cp =  NaN;
    ci =  NaN;
    k_c = NaN;
    td_c= NaN;
else
    qqi = zeros(len-1,1);
    for i = 2:len
        d1 = data(1:i-1,:);
        d2 = data(i:len,:);
        q1 = PosteriorHyper(d1,distr,a,knownp,[1 i-1]);
        q2 = PosteriorHyper(d2,distr,a,knownp,[i len]);
        qqi(i-1) = qqi(i-1) + MarginalModelLikelihood(d1,distr,q1,a,knownp,ops);
        qqi(i-1) = qqi(i-1) + MarginalModelLikelihood(d2,distr,q2,a,knownp,ops);
    end
    q = PosteriorHyper(data,distr,a,knownp,[1 len]);
    qqi = qqi - MarginalModelLikelihood(data,distr,q,a,knownp,ops);
    k_c = qqi;
    %==Calculate denominator==
    rtdur = times(len)-times(1);
    td_c = zeros(len-1,1);
    for i = 1:len-1
        td_c(i) = log(times(i+1)-times(i))-log(rtdur);
    end
    %==SBIC Model Complexity Correction==
    switch distr
        case 'linear'
            prms = (3.*size(a{5},2) + size(a{5},2).^2)./2;
        case 'multivariate normal'
            prms = (3.*ops + ops.^2)./2;
        case 'multiple linear'
            prms = ops.*(3.*size(a{5},2) + size(a{5},2).^2)./2;
        otherwise
            prms = sum(isnan(knownp));
    end
    if len > 2
        reltime = (0:len-1)'/(len-1);
        modif = zeros(len-1,1);
        dx = reltime(2);
        modif(1) = (dx.*(2-log(dx - dx.^2)) + log(1-dx)) - 2.*dx;
        dx = reltime(len) - reltime(len-1);
        modif(len-1) = (dx.*(2-log(dx - dx.^2)) + log(1-dx)) - 2.*dx;
        if len > 2
            for i = 2:len-2
                dx2 = reltime(i+1);
                dx1 = reltime(i);
                modif(i) = (dx2.*(2-log(dx2 - dx2.^2)) + log(1-dx2)) - (dx1.*(2-log(dx1 - dx1.^2)) + log(1-dx1)) - 2.*(dx2 - dx1); %Caculate finite integral
            end
        end
        modif = prms.*( modif.*len./2) ; %convert to SBIC correction
        td_c = td_c - modif;
    else
        modif = (0.5.*(2-log(0.5 - 0.5.^2)) + log(1-0.5));
        td_c = td_c - (modif).*(prms.*len);
    end
    %==Pick CP==
    k_c(isinf(k_c))=NaN;
    BF = exp(k_c + td_c).*weight;
    cp = find((k_c + td_c + log(weight))==max((k_c + td_c + log(weight))),1,'first');
    if isempty(cp) || cp==0
        cp = NaN;
        ci = NaN;
    end
    %==Calculate Incremental Confidence Interval==
    if SETTINGS.CONF
        ci = [cp cp];
        if sum(isinf(BF)) > 0
            BF(~isinf(BF)) = BF(~isinf(BF))./(2.*len);
            BF(isinf(BF)) = realmax./(2.*len);
        end
        BF(isnan(BF)) = 0;
        cb = (1-conf)./2;
        ct = 1-cb;
        BF = cumsum(BF./sum(BF));
        top = (BF>ct);
        bot = (BF<cb);
        if sum(bot)==0
            ci(1) = times(1) + (times(2)-times(1)).*(cb./BF(1));
		else
            q = find(bot,1,'last');
			if q < length(BF)
	            ci(1) = times(q) + (times(q+1)-times(q)).*((cb-BF(q))./(BF(q+1)-BF(q)));
			else
	            ci(1) = times(q) + (times(q+1)-times(q)).*((cb)./(BF(q)));				
			end
        end
        if sum(top)==0
            ci(2) = times(len-1) + (times(len)-times(len-1)).*(1-(cb./(1-BF(len-1))));
        else
            q = find(top,1,'first');
            if q > 1
                ci(2) = times(q) + (times(q+1)-times(q)).*((ct-BF(q-1))./(BF(q)-BF(q-1)));
            else
                ci(2) = times(1) + (times(2)-times(1)).*((ct)./(BF(1)));
            end
        end
    else
        ci = [NaN NaN];
     end
end
end
