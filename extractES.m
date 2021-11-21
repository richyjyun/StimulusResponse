%% Get evoked spikes with timestamp data

function [chns, codes, evokedtimes, evokeddelays, spiketimes] = extractES(spikes,stim)

[chns,codes,~] = findSpikeChnCode(spikes);
evokedtimes = {};
evokeddelays = {}; 
spiketimes = {};
for c = 1:length(chns)
      
    chn = chns(c); sc = codes(c);
    spike = spikes.timestamp(spikes.chan == chn & spikes.sortcode == sc);
    spike = spike+0.0005; % NEED TO ADD TO SHIFT FOR RIPPLE WINDOW
    
    % Must have at least 500 spikes
    minspikes = 500;
    if(length(spike) < minspikes)
        evokedtimes{c} = [];
        evokeddelays{c} = [];
        spiketimes{c} = spike;
        continue;
    end
    
    %% Remove spikes that are way too close to stimulation (<1.2ms)
    bins = [0,stim(1:end-1),max(stim(end),spike(end))+1];
    stimind = discretize(spike,bins);
    
    % time since last stim
    closest = bins(stimind)';
       
    delay = spike-closest;

    % time until next stim
    after = bins(stimind+1)';

    next = spike-after;
    
    artifact = 0.0012; 
            
    tooclose = delay < artifact | next > -artifact;
    
    %% Get PSTH
    rightcorbound = 0.8*mean(diff(stim));
    [cor,lags] = CrossCorr(stim, 'ts2',spike(~tooclose),'binsize', 0.0005,'lag',[-0.02,rightcorbound],'suppress_plot',1);

    %% Get threshold for evoked spikes
    peakbounds = [1.2,10]*0.001;
    
    tempi = find(lags>-0.002,1);
    premean = mean(cor(1:tempi));
    prestd = std(cor(1:tempi));
    
    thresh = premean+2*prestd;
    
    %% Find peak in PSTH and the bounds for evoked spikes
    [pks,loc] = findpeaks(cor);
    ind = find(lags(loc) > peakbounds(1) & lags(loc) < peakbounds(2));
    ind(cor(loc(ind)) < thresh) = [];
    if(isempty(ind)) % Check there are peaks at all
        evokedtimes{c} = [];
        evokeddelays{c} = [];
        spiketimes{c} = spike;
        continue;
    else
        pks = pks(ind); loc = loc(ind); % Limit peaks to bounds
        [maxpk,maxind] = max(pks); % Find maximum peak
        
        % Get bounds to find evoked spikes in
        evokedbounds = [-1.5,1.5]*0.001;
        shift = lags(loc(maxind));
        evokedbounds = evokedbounds + shift;
        if(evokedbounds(1) < 0)
            evokedbounds = [0,2]*0.001;
        end
    end
    
    bounds = evokedbounds;

    %% Find evoked spikes
    evoked = delay >= bounds(1) & delay <= bounds(2);
    
    if(sum(evoked) == 0)
        evokedtimes{c} = [];
        evokeddelays{c} = [];
        spiketimes{c} = spike;
        continue;
    end
    
    conddelay = delay(evoked);
    
    %% Save values
    spiketimes{c} = spike(~evoked & ~tooclose); 
    evokedtimes{c} = spike(evoked); 
    evokeddelays{c} = conddelay;
      
end

end