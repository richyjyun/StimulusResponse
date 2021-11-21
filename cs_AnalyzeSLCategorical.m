%% Analysis of changes over time 
% Categorical data, handled a bit differently

slname = 'SL.mat';
load(slname);
SL = SL(1:30);

% Load spanky data and format
spanky = 'C:\Users\Richy Yun\Dropbox\Fetz Lab\EvokedSpikes\Spanky Data\sorted_eSpikes_OverTime.mat';
load(spanky);
stimSpikes{5} = stimSpikes{3}{2};
stimSpikes{3} = stimSpikes{3}{1};
stimTimes{5} = stimTimes{3}{2};
stimTimes{3} = stimTimes{3}{1};

spankySL = struct([]);
for i = 3:length(stimSpikes)
    
    SD = struct();
    SD.spiketimes{1} = stimSpikes{i};
    SD.stim = stimTimes{i};
    SD.stimfreq = length(SD.stim)/(SD.stim(end)-SD.stim(1));
    
    % Evoked spikes
    [ES, dt, prc, norm] = getESpikes(stimSpikes{i}, stimTimes{i});
    SD.ES{1} = stimSpikes{i}(ES);
    SD.ES_delay{1} = dt;
    
    % Inhibition
    SD.IH = extractIH(SD.ES_delay,SD.spiketimes,SD.stim);
    
    % Evoked spike trends
    [SD.ESslope, SD.ESpval] = extractTrend(SD.ES,SD.stim,1);
    
    % Inhibition trends
    [SD.IHslope, SD.IHpval] = extractTrend(SD.IH,SD.stim,0);
    
    % Channel distance
    SD.dist = 0.4;
    
    if(isempty(spankySL))
        spankySL = SD;
    else
        spankySL(end+1) = SD;
    end
    
end


%% Compile all data
ESfreq = {};
ESchange = {};
IHfreq = {};
IHchange = {};

perlim = 0.03;

for s = 1:length(SL)
        
    temp = sign(SL(s).ESslope);
    ESper = cellfun(@(x) length(x)/length(SL(s).stim),SL(s).ES);
    bad = ESper < perlim | ESper > 0.95;
    temp(SL(s).ESpval >= 0.05) = 0;
    
    ESchange{s} = temp(~bad);
    ESfreq{s} = repmat(length(SL(s).stim)/(SL(s).stim(end)-SL(s).stim(1)),1,sum(~bad));
    
    temp = sign(SL(s).IHslope);
    bad = ESper < perlim | ESper > 0.95;
    temp(SL(s).IHpval >= 0.05) = 0;
    
    IHchange{s} = temp(~bad);
    IHfreq{s} = repmat(SL(s).stimfreq,1,sum(~bad));

end

allfreqs = cell2mat(ESfreq);
allESchange = cell2mat(ESchange);

low = allfreqs < 5;
high = allfreqs >= 15;
med1 = allfreqs >= 5 & allfreqs < 10;
med2 = ~low & ~high & ~med1;

Change{1} = [allESchange(low),1,1,0,0]; % adding lorde and spanky data
Change{2} = [allESchange(med1),-1,0,0];
Change{3} = allESchange(med2);
Change{4} = [allESchange(high),-1,1,-1,-1,-1,-1]; 

changepvals = [];
decreasepvals = [];
changepvalsf = [];
decreasepvalsf = [];
for i = 1:length(Change)
    for j = i+1:length(Change)
        n1 = sum(Change{i} == 0); n2 = sum(Change{i} ~= 0);
        n3 = sum(Change{j} == 0); n4 = sum(Change{j} ~= 0);
        
        [h,changepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,changepvals(i,j)] = crosstab(x1,x2);
        
        n1 = sum(Change{i} == -1); n2 = sum(Change{i} == 1);
        n3 = sum(Change{j} == -1); n4 = sum(Change{j} == 1);
        
        [h,decreasepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,decreasepvals(i,j)] = crosstab(x1,x2);
    end
end

ESchangeprob = [];
ESdecreaseprob = [];
for i = 1:length(Change)
    n1 = sum(Change{i} == 0); n2 = sum(Change{i} ~= 0);
    ESchangeprob(i) = n2/(n1+n2)*100;
    n1 = sum(Change{i} == -1); n2 = sum(Change{i} == 1);
    ESdecreaseprob(i) = n1/(n1+n2)*100;
end

%% Binned
bw = 4;
bins = 0:bw:max(allfreqs); bins = [bins,bins(end)+bw];
inds = discretize(allfreqs,bins);
changeprob = [];
decreaseprob = [];
for i = 1:(length(bins)-1)
    total = allESchange(inds==i);
    changeprob(i) = sum(total~=0)/length(total);
    decreaseprob(i) = sum(total==-1)/sum(total~=0);
end

figure; plot(bins(1:end-1),changeprob);
hold on; plot(bins(1:end-1),decreaseprob);

figure; scatter(bins(1:end-1),changeprob);
hold on; scatter(bins(1:end-1),decreaseprob);

%% Inhibition
allfreqs = cell2mat(IHfreq);
allESchange = cell2mat(IHchange);

low = allfreqs < 5;
high = allfreqs >= 15;
med1 = allfreqs >= 5 & allfreqs < 10;
med2 = ~low & ~high & ~med1;

Change{1} = [allESchange(low),1,-1,0,0];
Change{2} = [allESchange(med1),0,0,0];
Change{3} = allESchange(med2);
Change{4} = [allESchange(high),0,1,0,0,0,0]; 

changepvals = [];
decreasepvals = [];
changepvalsf = [];
decreasepvalsf = [];
for i = 1:length(Change)
    for j = i+1:length(Change)
        n1 = sum(Change{i} == 0); n2 = sum(Change{i} ~= 0);
        n3 = sum(Change{j} == 0); n4 = sum(Change{j} ~= 0);
        
        [h,changepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,changepvals(i,j)] = crosstab(x1,x2);
        
        n1 = sum(Change{i} == -1); n2 = sum(Change{i} == 1);
        n3 = sum(Change{j} == -1); n4 = sum(Change{j} == 1);
        
        [h,decreasepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,decreasepvals(i,j)] = crosstab(x1,x2);
    end
end

IHchangeprob = [];
IHdecreaseprob = [];
for i = 1:length(Change)
    n1 = sum(Change{i} == 0); n2 = sum(Change{i} ~= 0);
    IHchangeprob(i) = n2/(n1+n2)*100;
    n1 = sum(Change{i} == -1); n2 = sum(Change{i} == 1);
    IHdecreaseprob(i) = n1/(n1+n2)*100;
end

%% Figure
set(0,'defaultAxesFontSize',12)

figure; subplot(2,2,1); plot(ESchangeprob); hold on; plot(ESdecreaseprob)
xticks([1:4]); xlim([0.5,4.5]); xlabel('Stim Frequency')
xticklabels({'0-5Hz','5-10Hz','10-15Hz','15-20Hz'});
ylabel('Probability (%)'); title('Changes in Evoked Spikes');
legend({'Change','Decrease'},'box','off');  box off;

subplot(2,2,2);
plot(IHchangeprob); hold on; plot(IHdecreaseprob)
xticks([1:4]); xlim([0.5,4.5]); xlabel('Stim Frequency')
xticklabels({'0-5Hz','5-10Hz','10-15Hz','15-20Hz'});
title('Changes in Inhibition'); box off;













