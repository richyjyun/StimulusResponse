%% Get the inhibitory response
% Basically finding the next spike after stimulation that isn't an evoked
% spike

function inhib = extractIH(delays,spikes,stim)
    
inhib = cell(1,length(spikes));
for c = 1:length(spikes)
    
    if(isempty(delays{c}))
        mintime = 0;
    else
        mintime = mean(delays{c});
    end
    
    temp = discretize(spikes{c},stim);
    nextspike = nan(1,length(stim));
    for s = 1:(length(stim)-1)
        temp = find(spikes{c}>(stim(s)+mintime) & spikes{c}<(stim(s+1)),1);
        if(~isempty(temp))
            nextspike(s) = temp;
        end
    end
    
    good = ~isnan(nextspike);
    
    temp1 = spikes{c}(nextspike(good));
    temp2 = stim(good);
    nextspike(good) =  temp1(:) - temp2(:);
    
    inhib{c} = nextspike;
    
end

end
