slname = 'R:\Yun\Kronk\Ripple\SL_ESFix.mat'; 
load(slname);

%% Consolidate channels and spikes
Chns = []; 
Codes = [];
Stimchn = [];
for s = 1:length(SL)
    Stimchn(s) = SL(s).stimchn;
    chns = SL(s).chns;
    codes = SL(s).codes;
    for c = 1:length(chns)
       inds = find(Chns==chns(c));
       if isempty(inds) || ~any(Codes(inds)==codes(c))
           Chns(end+1) = chns(c);
           Codes(end+1) = codes(c);
       end
    end
end

% 148 unique spikes across 69 channels
% 16 unqiue stimulated channels

%% Tagging each spike in each session to track between sessions
for s = 1:length(SL)
    chns = SL(s).chns;
    codes = SL(s).codes;
    tag = [];
    for c = 1:length(chns)
        tag(c) = find(Chns==chns(c) & Codes==codes(c));
    end
    SL(s).tag = tag;
    
end

%% Consistency of decrease / increase per spike
ESChange = cell(1,length(Chns));
IHChange = cell(1,length(Chns));
Freq = cell(1,length(Chns));
Stim = cell(1,length(Chns));
for s = 1:length(SL)
    tag = SL(s).tag;
    for t = 1:length(tag)
        eschange = 0;
        if SL(s).ESpval(t) < 0.05
            eschange = sign(SL(s).ESslope(t));
        end
        ESChange{tag(t)} = [ESChange{tag(t)},eschange];
        
        ihchange = 0;
        if SL(s).IHpval(t) < 0.05
            ihchange = sign(SL(s).IHslope(t));
        end
        IHChange{tag(t)} = [IHChange{tag(t)},ihchange];
    
        Freq{tag(t)} = [Freq{tag(t)},SL(s).stimfreq];
        Stim{tag(t)} = [Stim{tag(t)},SL(s).stimchn];

    end
    
end

i = 1;
jitter = 2*rand(1,length(Freq{i}))-1;
figure; scatter(Freq{i}+jitter,ESChange{i});

i = 4;
jitter = 2*rand(1,length(Freq{i}))-1;
figure; scatter(Freq{i}+jitter,IHChange{i});

%% Determine if position in the array affects changes

allChns = unique(Chns);
grid = nan(10,10);
for i = 1:length(allChns)
    ind = find(Chns==allChns(i));
    temp = ESChange(ind);
    temp = cell2mat(temp);
    
    [c, r] = GetKronkChannelPosition(allChns(i));
    
    val = mean(temp);%sum(temp==1)./sum(temp==1|temp==-1);
    if(sum(temp==1|temp==-1) == 0)
        val = 0;
    end
    grid(c,r) = val;
end
figure; imagesc(grid,'AlphaData',~isnan(grid)); colorbar;

% Grid
hold on;
edges = 0.5:1:10.5;
for i = 1:11
    if(i==1 || i==11)
        plot([edges(2),edges(end-1)],[edges(i),edges(i)],'k');
    else
        plot([edges(1),edges(end)],[edges(i),edges(i)],'k');
    end
end
for i = 1:11
    if(i==1 || i==11)
        plot([edges(i),edges(i)],[edges(2),edges(end-1)],'k');
    else
        plot([edges(i),edges(i)],[edges(1),edges(end)],'k');
    end
end

%% Determine if position relative to the stimulus channel affects changes
grid = cell(19,19);
for i = 1:length(ESChange)
    temp = ESChange{i};
    sc = []; sr = [];
    for s = 1:length(Stim{i})
        [sc(s), sr(s)] = GetKronkChannelPosition(Stim{i}(s));
    end
    [c,r] = GetKronkChannelPosition(Chns(i));
    
    col = sc-c + 10; row = sr-r + 10;
    
    for s = 1:length(temp)
        grid{row(s),col(s)} = [grid{row(s),col(s)}, temp(s)];
    end
end

temp = cellfun(@(x) sum(x==1)./sum(x==-1 | x==1),grid);
temp2 = cellfun(@(x) sum(x==-1|x==1),grid);

figure; imagesc(temp,'AlphaData',~isnan(temp)); colormap jet;
colorbar;
hold on; fill([9.5,10.5,10.5,9.5],[9.5,9.5,10.5,10.5],'k');


%% Locations of all stim channels
temp = cell2mat(Stim);
temp = unique(temp);
grid = nan(10,10);
for s = 1:length(temp)
    [c,r] = GetKronkChannelPosition(temp(s));
    grid(r,c) = 0;
end

figure; imagesc(grid,'AlphaData',~isnan(grid));


%% Does each individual stimulus have a more likely chance of evoking a population of spikes?
% Comparing counts
tempSL = SL(28);
counts = zeros(1,length(tempSL.stim));
counts_c = zeros(1,length(tempSL.stim));
sum(cellfun(@(x) ~isempty(x), tempSL.ES))
for i = 1:length(tempSL.ES)
    if(~isempty(tempSL.ES{i}))
        ind = discretize(tempSL.ES{i},tempSL.stim);
        counts(ind) = counts(ind)+1;
        
        temp = randperm(length(counts));
        temp = temp(1:length(ind));
        counts_c(temp) = counts_c(temp)+1;
        
    end   
end

figure; subplot(2,1,2);
histogram(counts,'edgealpha',0);
ylabel('Stim Count'); xlabel('Simultaneous Evoked Spikes')
hold on; histogram(counts_c,'edgealpha',0);
[~, p] = kstest2(counts, counts_c); box off;
set(gca, 'fontsize', 10);

%% Pairs of spikes with cross correlation w.r.t. distance
corr = {};
R = [];
P = [];
d = [];
m = [];
sim = [];
p_sim = [];
for s = 1:length(SL)
    tempSL = SL(s);
    for i = 1:(length(tempSL.ES)-1)
        ES1 = zeros(1,length(tempSL.stim));
        if(~isempty(tempSL.ES{i}))
            ind = discretize(tempSL.ES{i},tempSL.stim);
            ES1(ind) = 1;
        else
            continue;
        end
        for j = (i+1):length(tempSL.ES)
            ES2 = zeros(1,length(tempSL.stim));
            if(~isempty(tempSL.ES{i}))
                ind = discretize(tempSL.ES{j},tempSL.stim);
                ES2(ind) = 1;
                [temp,lags] = xcorr(ES1,ES2,10);
                if(all(temp==0))
                    continue;
                end
                corr{end+1} = temp;
                
                [r,p] = corrcoef(ES1,ES2);
                R(end+1) = r(1,2);
                P(end+1) = p(1,2);
                
                sim(end+1) = jaccard(ES1,ES2);
                
                [c1,r1] = GetKronkChannelPosition(tempSL.chns(i));
                [c2,r2] = GetKronkChannelPosition(tempSL.chns(j));
                
                d(end+1) = pdist2([c1,r1],[c2,r2])*0.4;
                m(end+1) = abs(c1-c2)+abs(r1-r2)*0.4;
                
            else
                continue;
            end
            
        end
    end
end

allcorr = cell2mat(corr')';
% figure; plot(allcorr);
% figure; plot(zscore(allcorr));

bad = d==0; 
jw = 0.1; 
jitter = jw*rand(1,length(d))-(jw/2);
sig = P<0.01;
subplot(2,2,1);
scatter(d(~bad&~sig)+jitter(~bad&~sig),R(~bad&~sig),'k.');
hold on; scatter(d(~bad&sig)+jitter(~bad&sig),R(~bad&sig),'r.');
ylabel('R'); xlabel('Distance (mm)');
xlim([0,5]); xl = xlim; hold on; plot(xl,[0,0],'k--');
set(gca, 'fontsize', 10);

% Binned % of covarying vs not
bins = 0.25:0.5:4;
inds = discretize(d(~bad), bins);
temp = R(~bad);

ratio = [];
for i = 1:max(inds)
    ratio(end+1) = sum(sig(inds==i)) / length(sig(inds==i)) * 100;
end
p = [];
for i = 1:max(inds)-1
    for j = i:max(inds)
        
        n1 = sum(sig(inds==i)); n2 = sum(~sig(inds==i));
        n3 = sum(sig(inds==j)); n4 = sum(~sig(inds==j));
        
        [~,p(i,j)] = fishertest([n1,n2; n3,n4]);
    end
end

subplot(2,2,2); 
plot([0.5:0.5:3.5], ratio, 'k');
hold on; scatter([0.5:0.5:3.5], ratio, 'k', 'filled');
xlim([0.25, 3.75]); box off;
ylim([10,25])
xlabel('Distance (mm)');
ylabel('Percent');
set(gca, 'fontsize', 10);











