%% Linear and exponential regression to see if changes over time are significant

function [slope,pval] = extractTrend(data,stimtimes,ES)

slope = zeros(1,length(data));
pval = nan(1,length(data));
    
% setting up bins
bw = 1; bins = stimtimes(1):bw:stimtimes(end);
Time = bins(1:end-1)+bw/2;

% function for exponential fit
modelfun = @(b,x) b(1) + b(2).*(1-b(3)).^x(:,1);

% accepted p-value
pvallim = 0.05;

% for debugging
plim = [];
pvarlim = [];

warning ('off','all');
for c = 1:length(data)
    if(isempty(data{c}))
        slope(c) = nan;
        continue;
    end
    
    %% Calculate binned data over time
    
    % For evoked spikes
    if ES
        evokedcount = histcounts(data{c},bins);
        stimcount = histcounts(stimtimes,bins);
        binned = evokedcount./stimcount;
        
        % For inhibition
    else
        stimind = discretize(stimtimes,bins);
        binned = arrayfun(@(x) nanmedian((data{c}(stimind==x))) ,1:max(stimind));
    end

    % For values that are still nan interpolate
    bad = isinf(binned);
    binned(bad) = nan; time = Time; 
    
    % There needs to be enough ES to properly do regression
    if(sum(binned==0 | isnan(binned)) > length(binned)/4)
        continue;
    end
    
    binned = smooth(binned,30); % smooth over 10 seconds, so essentially 10 second bins with 1 second steps
    
    %% Regression
    % Linear regression
    try
        mdl = fitlm(time',binned','RobustOpts','huber');
    catch
        continue
    end
    f = [mdl.Coefficients{2,1},mdl.Coefficients{1,1}];
    
    % Fitted line
    y = polyval(f,time);
    
    p = mdl.Coefficients{2,4};
    
    if p<pvallim
        plim = [plim,c];
        if(abs((y(end)-y(1)))>std(binned))
            pvarlim = [pvarlim,c];
            slope(c) = f(1);
            pval(c) = p;
        end
    else
        
        % Exponential regression (decrease)
        try
            tbl = table(time',binned');
            beta0 = [mean(binned(1:5)), 2, 0.01];
            mdl = fitnlm(tbl,modelfun,beta0);
            
            % Fitted curve
            y = mdl.Fitted;
            
            % If exponent slopes the wrong way, keep going
            if mdl.Coefficients{2,1}<0 || mdl.Coefficients{3,1}>1
                continue;
            end
            
            p1 = mdl.Coefficients{1,4};
            p2 = mdl.Coefficients{2,4};
            p3 = mdl.Coefficients{3,4};
            
            pExp = max([p1,p2,p3]);
            if(all(~isnan([p1,p2,p3])) && all([p1,p2,p3]~=0) && pExp<pvallim)
                if(abs((y(end)-y(1)))>std(binned))
                    slope(c) = f(1);
                    pval(c) = pExp;
                end
            end
            
        catch
        end
        
    end
    
end

warning ('on','all');


end