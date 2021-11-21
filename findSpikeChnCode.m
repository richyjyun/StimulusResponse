%% Finds all channels and sortcodes that were saved for single units
% Takes sorting snip struct from TDT as input, returns channels and codes

function [chns,codes,counts] = findSpikeChnCode(Snips)

channels = unique(Snips.chan);
chns = []; codes = []; counts = [];
for c = 1:length(channels)
    sc = unique(Snips.sortcode(Snips.chan==channels(c)));
    sc(sc<1 | sc>5) = [];
    if(~isempty(sc))
        codes = [codes;sc];
        chns = [chns;double(channels(c)).*ones(length(sc),1)]; 
        for i = 1:length(sc)
            counts = [counts;sum(Snips.chan==channels(c) & Snips.sortcode == sc(i))];
        end
    end
end

end