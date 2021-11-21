%% Find covariance of evoked spikes - i.e. do different units get evoked by the same stimuli

function Cov = extractCov(es,stim)
 
evokedstim = cell(1,length(es));
for i = 1:length(es)
    temp = discretize(es{i}, stim);
    temp(isnan(temp)) = [];
    evokedstim{i} = temp;
end

Cov = zeros(length(es),length(es));
for i = 1:(length(es)-1)
    if(isempty(evokedstim{i}))
        continue;
    end
    for j = (i+1):length(es)
        if(isempty(evokedstim{j}))
            continue;
        end
        temp1 = zeros(1,length(stim));
        temp2 = zeros(1,length(stim));
        temp1(evokedstim{i}) = 1;
        temp2(evokedstim{j}) = 1;
        [cor,lags] = xcorr(temp1,temp2, 20, 'coeff');
        [Corr, ind] = max(cor);
        
        % Make sure maximum happenes at zero, the same stim
        if(lags(ind)==0)
            tempcor = cor; tempcor(ind-1:ind+1) = [];
            % The correlation at zero must be large relative to the rest
            if(Corr > mean(tempcor)+2.5*std(tempcor))
                Cov(i,j) = 1;
            else
                Cov(i,j) = -1;
            end
        else
            Cov(i,j) = -1;
        end
    end
end

% Reflect across diagonal to make it easier later
Cov = Cov+Cov';

end