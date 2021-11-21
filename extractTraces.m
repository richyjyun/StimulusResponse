%% Extracting spike traces to save for cell type analysis 

function traces = extractTraces(path,filename,rawspikes,spiketimes)       

[chns,codes,~] = findSpikeChnCode(rawspikes);
prevchn = 0;

traces = cell(1,length(chns));
for c = 1:length(chns)
    
    if(chns(c) ~= prevchn)
        t = LoadSpikeTrace(path,filename,chns(c));
        prevchn = chns(c);
        if(isempty(t) || all(t.sortcode==0)), continue; end
    end
    
    spike = find(t.sortcode == codes(c));
    
    [~,keep,~] = intersect(t.timestamp(spike)+0.0005, spiketimes{c});
    
    traces{c} = mean(t.snips(:,spike(keep)),2);
    
end

end