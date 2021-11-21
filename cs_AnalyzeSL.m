%% Broad analyses of stimulus evoked spikes and the inhibitory response

slname = 'SL.mat';
load(slname);
SL = SL(1:30);

%% Evoked spike example
s = 46; c = 68;
set(0,'defaultAxesFontSize',12)

figure; 
spikes = sort([SL(s).spiketimes{c};SL(s).ES{c}]);
[cor,lags] = CrossCorr2_Clean(SL(s).stim,spikes+1e-3,300,2000,0.05,1,2);
zerolag = find(lags==0); 
PSTH = smoothdata(cor,'gaussian',10);
PSTH = (cor+PSTH)/2; PSTH(zerolag:zerolag+1) = 0;
basewin = [find(lags==-0.02), find(lags==-0.002)];

bar(lags*1000,PSTH,1,'k');

base = PSTH(basewin(1):basewin(2));
upperthresh = mean(base)+2*std(base);
lowerthresh = mean(base)-2*std(base);

leg = [];

ylim([0,0.15]);

xl = xlim;
hold on; 
leg(1) = plot(xl,[upperthresh,upperthresh],'b--','linewidth',1.5);
leg(2) = plot(xl,[lowerthresh,lowerthresh],'r--','linewidth',1.5);

eslim = [1.25,5.25];
yl = ylim;
leg(3) = patch([eslim(1),eslim,fliplr(eslim)],[yl,fliplr(yl),yl(1)],...
    'b','facealpha',0.2,'edgealpha',0);
xlim([-20,20]);

hold on; plot([0,0],yl,'k--');

legend(leg,{'Upper Threshold','Lower Threshold','Evoked Spikes'},'box','off');

xlabel('Time since stim (ms)');
ylabel('Probability');

box off;

%% Inhib Example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SD = SL(27); 
c = 49; 
spike = sort([SD.ES{c};SD.spiketimes{c}]);
raster = {};
figure; 
[cor,lags] = CrossCorr(SD.stim, 'ts2',spike,'binsize', 0.001,'lag',[-0.02,0.16],'suppress_plot',1);
set(0,'defaultAxesFontSize',10)
subplot(6,1,[1:2]); cor(lags>0 & lags<0.002) = 0;
bar(lags,cor,1,'k','edgealpha',0,'facealpha',1); box off;
set(gca,'FontSize',10)
yl = ylim; hold on; plot([0,0],yl,'k--','linewidth',2);
xlim([-0.02,0.16]); ylabel('Probability'); xticks([]);
subplot(6,1,[3:6]);
n = 1500;
h = [];
for i = 1:n
    temp = spike-SD.stim(i);
    temp(temp<-0.02 | temp>0.16) = [];
    temp(temp>0 & temp<0.002) = [];
    raster{i} = temp;
    hold on; h(1) = scatter(temp,(SD.stim(i)-SD.stim(1))*ones(1,length(temp)),'.','markeredgecolor',[0.7,0.7,0.7]);
end
hold on; h(3) = scatter(SD.IH{c}(1:n),SD.stim(1:n)-SD.stim(1),'k.');
estim = discretize(SD.ES{c},SD.stim);
hold on; h(2) = scatter(SD.ES_delay{c},SD.stim(estim)-SD.stim(1),'b.');
bw = 10; bins = SD.stim(1):bw:(SD.stim(n)+bw);
ind = discretize(SD.stim(1:n),bins);
temp = SD.IH{c}(1:n);
avg = arrayfun(@(x) nanmedian(temp(ind==x)), 1:max(ind));

x = SD.stim(1):SD.stim(n);
vq = interp1(bins(1:end-1)+bw/2,avg,x,'pchip');

hold on; h(4) = plot(vq,x-SD.stim(1),'r','linewidth',2.5);
xlim([-0.02,0.16]); ylim([0,SD.stim(n)-SD.stim(1)]);
yl = ylim; hold on; plot([0,0],yl,'k--','linewidth',2);
xlabel('Time since stim (s)');
ylabel('Experiment time (s)');
set(gca, 'YDir', 'reverse');
set(gca,'FontSize',10)


%% Distance figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distance vs % evoked
set(0,'defaultAxesFontSize',12)

D = [];
ES = [];
for s = 1:length(SL)
    ind = cellfun(@(x) ~isempty(x), SL(s).ES);
    temp = SL(s).ES(ind);
    prob = cellfun(@(x) length(x)/length(SL(s).stim),temp);
    ES = [ES,prob];
    D = [D,SL(s).dist(ind)];
end

bins = 0.25:0.25:3;
ind = discretize(D,bins);
avg = arrayfun(@(x) nanmedian(ES(ind==x)*100), 1:max(ind));

figure; s1 = subplot(2,2,1)
plot(avg,'color',[0.4,0.4,0.4]); hold on;
boxplot(ES*100,ind,'notch','on','symbol','w'); ylim([-2,85]);
set(findobj(gca,'type','line'),'linew',1.5); 
ylabel('% Evoked'); xlabel('Distance (mm)')
xticks(1.5:2:length(bins));
xticklabels(bins(2:2:length(bins))); xtickangle(45);
box off;

% Inhibition strength vs distance
D = [];
IH = [];
for s = 1:length(SL)
    ind = cellfun(@(x) ~isempty(x), SL(s).ES);
    temp = cellfun(@(x) nanmean(x),SL(s).IH(ind));
    IH = [IH,temp];
    D = [D,SL(s).dist(ind)];
end

bins = 0.25:0.25:3;
ind = discretize(D,bins);
avg = arrayfun(@(x) nanmedian(IH(ind==x)*1000), 1:max(ind));
s2 = subplot(2,2,2); plot(avg,'color',[0.4,0.4,0.4]); hold on;
boxplot(IH*1000,ind,'notch','on','symbol','w'); ylim([0,180]); 
set(findobj(gca,'type','line'),'linew',1.5); 
ylabel('Inhibition (ms)'); xlabel('Distance (mm)')
xticks(1.5:2:length(bins));
xticklabels(bins(2:2:length(bins))); xtickangle(45);
box off;
s2.Position(2:4) = s1.Position(2:4)

% Evoked spike timing 
D = [];
ES = [];
for s = 1:length(SL)
    ind = cellfun(@(x) ~isempty(x), SL(s).ES);
    temp = SL(s).ES(ind);
    for i = 1:length(temp)
        bin = discretize(temp{i},SL(s).stim);
        good = ~isnan(bin);
        ES(end+1) = nanmean(temp{i}(good)'-SL(s).stim(bin(good)));
    end
    D = [D,SL(s).dist(ind)];
end

subplot(2,2,[3,4]);

bins = 1:0.5:20; bins = bins; lim = 1;
hold on; histogram(ES(D<=lim)*1000,bins,'normalization','probability',...
    'edgealpha',0,'facealpha',0.7);
hold on; histogram(ES(D>lim)*1000,bins,'normalization','probability',...
    'edgealpha',0,'facealpha',0.7);

c = get(gca,'colororder');
counts = histcounts(ES(D<lim)*1000,bins); counts = smooth(counts./sum(counts),3);
x = bins(1:end-1)+0.25; xq = (1:0.01:20); vq = interp1(x,counts,xq,'pchip');
hold on; 
leg(1) = plot(xq,vq,'color',[c(1,:),0.5],'linewidth',3);
counts = histcounts(ES(D>=lim)*1000,bins); counts = smooth(counts./sum(counts),3);
x = bins(1:end-1)+0.25; xq = (1:0.01:20); vq = interp1(x,counts,xq,'pchip');
hold on; 
leg(2) = plot(xq,vq,'color',[c(2,:),0.5],'linewidth',3);
yl = ylim; ylim([0,max(yl)]);

legend(leg,{'<=1mm','>1mm'},'box','off'); 

xlim([1,20]);
xlabel('Evoked spike timing (ms)');
ylabel('Percent')

ytix = get(gca, 'YTick')
set(gca, 'YTick',ytix, 'YTickLabel',ytix*100)

%% Change with respect to frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultAxesFontSize',10)

% All plots below seem good. 
% Evoked spikes
ES = {};
ESslope = [];
ESper = [];
EStime = [];
ESpval = [];
perlim = 0.05;
for i = 1:length(SL)
    ES{i} = sign(SL(i).ESslope);
    

    ESpval = [ESpval,SL(i).ESpval];
    
    temp = cellfun(@(x) length(x)/length(SL(i).stim),SL(i).ES);
    bad = temp<perlim;
    temp = temp(~bad);
    ESper = [ESper,temp];
    
        ESslope = [ESslope,SL(i).ESslope(~bad)];
    
    temp = cellfun(@nanmedian, SL(i).ES_delay);
    EStime = [EStime,temp];
    
    ES{i} = ES{i}(~bad);
end

freq = extractfield(SL,'stimfreq');
groups = 2*ones(1,length(freq));
groups(freq<=7) = 1; groups(freq>14) = 3;

% Change vs no change
change = cellfun(@(x) sum(~isnan(x) & ~(x==0))/sum(~isnan(x)), ES);

figure; subplot(1,2,1);
boxplot(change*100,groups,'notch','on','symbol','w');
% hold on; xl = xlim; plot(xl,[50,50],'k--');
xticklabels({'Low','Medium','High'});
ylabel('Percentage of spikes'); xlabel('Stimulation Frequency');
title({'Evoked Spikes','Change over time'})
ylim([0,80]); set(gca,'box','off');
% sigstar({[1,3]},[.023],0,20,0);
set(findobj(gca,'type','line'),'linew',2);  ylim([0,100])

% Decrease probability
decrease = cellfun(@(x) sum(x==-1)/sum(~isnan(x) & x~=0), ES);

subplot(1,2,2);
boxplot(decrease*100,groups,'notch','on','symbol','w');
xticklabels({'Low','Medium','High'});
ylabel('Percentage of spikes'); xlabel('Stimulation Frequency');
title({'Evoked Spikes','Probability of decreasing'})
ylim([0,80]); set(gca,'box','off');
text(3,10,'**','horizontalalignment','center','fontsize',20)
set(findobj(gca,'type','line'),'linew',2);  ylim([0,100])
hold on; plot([.5,3.5],[50,50],'k--'); 

% Distance
D = {};
S = {};
for s = 1:length(SL)
    D{s} = SL(s).dist(~isnan(SL(s).ESslope));
    S{s} = sign(SL(s).ESslope(~isnan(SL(s).ESslope)));
end

bin = 0:0.5:3.5;
figure;
matches = [0,1,-1];
titles = {'Change','Increase','Decrease'};
for m = 1:3
    
    allavg = [];
    for i = 1:length(D)
        ind = discretize(D{i},bin);
        if(m == 1)
            avg = arrayfun(@(x) sum(S{i}(ind==x)~=matches(m))/sum(ind==x), unique(ind));
        else
            avg = arrayfun(@(x) sum(S{i}(ind==x)==matches(m))/sum(ind==x), unique(ind));
        end
        n = min(length(avg),length(bin)-1);
        allavg(i,1:n) = avg(1:n);
    end
    subplot(1,3,m); 
    plot(nanmedian(allavg),'color',[0.5,0.5,0.5]); hold on;
    boxplot(allavg,'symbol','w','notch','on');
    xticklabels(bin(1:end-1)+0.25); xlabel('Distance (mm)');
    ylabel('Probability');
    set(findobj(gca,'type','line'),'linew',1.5)
    ylim([-.1,1.05])
    xtickangle(45);
    title(titles{m});
end

% Evoked probability vs change
bad = isnan(ESslope);
ESslope(bad) = [];
ESpertemp = ESper; ESpertemp(bad) = [];
EStimetemp = EStime; EStimetemp(bad) = [];

groups = 1*(sign(ESslope)==-1) + 2*(sign(ESslope)==1) + 3*(ESslope==0);

figure; subplot(1,2,1);
boxplot(ESpertemp,groups);

subplot(1,2,2);
boxplot(EStimetemp,groups);

% Inhibition
IH = {};
for i = 1:length(SL)
    IH{i} = sign(SL(i).IHslope);
end
freq = extractfield(SL,'stimfreq');
groups = 2*ones(1,length(freq));
groups(freq<=7) = 1; groups(freq>14) = 3;

% Change vs no change
change = cellfun(@(x) sum(~isnan(x) & ~(x==0))/sum(~isnan(x)), IH);

figure; subplot(1,2,1);
boxplot(change*100,groups,'notch','on','symbol','w');
% hold on; xl = xlim; plot(xl,[50,50],'k--');
xticklabels({'Low','Medium','High'});
ylabel('Percentage of spikes'); xlabel('Stimulation Frequency');
title({'Inhibition','Change over time'})
ylim([0,80]); set(gca,'box','off');
set(findobj(gca,'type','line'),'linew',2);  ylim([0,100])
sigstar({[2,3],[1,3]},[.001,.02],0,20,0);

% Decrease probability
decrease = cellfun(@(x) sum(x==-1)/sum(~isnan(x) & x~=0), IH);

subplot(1,2,2);
boxplot(decrease*100,groups,'notch','on','symbol','w');
xticklabels({'Low','Medium','High'});
ylabel('Percentage of spikes'); xlabel('Stimulation Frequency');
title({'Inhibition','Probability of decreasing'})
ylim([0,80]); set(gca,'box','off');
text(2,10,'*','horizontalalignment','center','fontsize',20)
set(findobj(gca,'type','line'),'linew',2);  ylim([0,100])
hold on; plot([.5,3.5],[50,50],'k--'); 

% Distance
D = {};
S = {};
for s = 1:length(SL)
    D{s} = SL(s).dist(~isnan(SL(s).IHslope));
    S{s} = sign(SL(s).IHslope(~isnan(SL(s).IHslope)));
end

bin = 0:0.5:3.5;
figure;
matches = [0,-1,1];
titles = {'Change','Decrease','Increase'};
for m = 1:3
    
    allavg = [];
    for i = 1:length(D)
        ind = discretize(D{i},bin);
        if(m == 1)
            avg = arrayfun(@(x) sum(S{i}(ind==x)~=matches(m))/sum(ind==x), unique(ind));
        else
            avg = arrayfun(@(x) sum(S{i}(ind==x)==matches(m))/sum(ind==x), unique(ind));
        end
        n = min(length(avg),length(bin)-1);
        allavg(i,1:n) = avg(1:n);
    end
    subplot(1,3,m); 
    plot(nanmedian(allavg),'color',[0.5,0.5,0.5]); hold on;
    boxplot(allavg,'outliersize',0.001,'notch','on');
    xticklabels(bin(1:end-1)+0.25); xlabel('Distance (mm)');
    ylabel('Probability');
    set(findobj(gca,'type','line'),'linew',1.5)
    ylim([-.1,1.05])
    xtickangle(45);
    title(titles{m});
end

% Evoked probability vs change
bad = isnan(IHslope);
IHslope(bad) = [];
ESpertemp = ESper; ESpertemp(bad) = [];
EStimetemp = EStime; EStimetemp(bad) = [];

groups = 1*(sign(IHslope)==-1) + 2*(sign(IHslope)==1) + 3*(IHslope==0);

figure; subplot(1,2,1);
boxplot(ESpertemp,groups);

subplot(1,2,2);
boxplot(EStimetemp,groups);

%% Spike type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of spike width
% RMS or integration seems to give best distribution. Just calculating
% width is too small resolution.
% PCA doesn't do much
widths = [];
times = [];
ISI = [];
ES = [];
for s = 1:length(SL)
   traces = SL(s).traces;
   temp = cellfun(@(x) find(x==max(x),1)-find(x==min(x),1), traces,'uniformoutput',false);
   temp(cellfun(@isempty,temp)) = []; 
   times = [times,cell2mat(temp)];
   
   temp = cellfun(@(x) sum(abs(x))./(max(abs(x))),traces,'uniformoutput',false);
   temp(cellfun(@isempty,temp)) = [];
   widths = [widths,cell2mat(temp)];

   temp = cellfun(@(x) nanmedian(diff(x)), SL(s).spiketimes);
   ISI = [ISI,temp];
   
   es = cellfun(@(x) length(x)/length(SL(s).stim), SL(s).ES);
   ES = [ES,es];

end
bad = times<0 | times/30>1;
traces = SL(1).traces;
set(0,'defaultAxesFontSize',12)
figure; 
subplot(5,2,1); plot((0:51)/30,traces{21}*2,'linewidth',2); 
xlabel('ms'); yticks([]); box off; xlim([0,51/30]);
hold on; plot((0:51)/30,traces{1},'linewidth',2);  
xlabel('ms'); yticks([]); box off; xlim([0,51/30]);
legend({'RS','FS'})

subplot(5,2,2);  histogram(widths(~bad),100,'normalization',...
    'probability','edgealpha',0,'facecolor','k','facealpha',1);
ylabel('Probability'); xlabel('Spike area under curve (a.u.)');
yl = ylim; hold on; plot([12,12],yl,'r--'); box off;


[coeff,score,latent,tsquared,explained,mu] = pca(shapes');%,'Algorithm','eig');
Vr = score;
figure; scatter3(Vr(widths<7.5,1),Vr(widths<7.5,2),Vr(widths<7.5,3),'.');
hold on; scatter3(Vr(widths>=7.5,1),Vr(widths>=7.5,2),Vr(widths>=7.5,3),'.');

% PCA
widths = [];
for s = 1:length(SL)
   temp = cell2mat(SL(s).traces);
   temp = temp./(max(temp)-min(temp));
   
   widths = [widths,temp];
end
[coeff,score,latent,tsquared,explained,mu] = pca(widths');%,'Algorithm','eig');
variance = cumsum(explained/sum(explained));
r = find(variance > 0.90,1);
Vr = score(:,1:r);

idx = kmeans(Vr,2,'Replicates',50,'MaxIter',1000);
figure;
for i = 1:2
    hold on; scatter3(Vr(idx==i,1),Vr(idx==i,2),Vr(idx==i,3),'.');
end

% ISI vs evoked prob
bad = isnan(ESslope) | isempty(ESslope) | isnan(ISI);
EStemp = ESslope(~bad);
ISItemp = ISI(~bad);

figure; boxplot(EStemp,ISItemp>0.08);


% Width vs evoked probability
allwidths = [];
ES = []; IH = [];
dist = [];
ESchange = []; IHchange = [];
for s = 1:length(SL)
    temp = SL(s).traces;
    temp = cellfun(@(x) sum(abs(x))./max(abs(x)),temp,'uniformoutput',false);
    bad = cellfun(@isempty,temp);
    temp(bad) = [];
    allwidths = [allwidths,cell2mat(temp)];
    
    es = cellfun(@(x) length(x)/length(SL(s).stim), SL(s).ES);
    es(bad) = [];
    ES = [ES,es];
    
    ih = cellfun(@(x) nanmedian(x), SL(s).IH);
    ih(bad) = [];
    IH = [IH,ih];
    
    d = SL(s).dist;
    d(bad) = [];
    dist = [dist,d];
    
    temp = sign(SL(s).ESslope);
    ESchange = [ESchange,temp];
    
    temp = sign(SL(s).IHslope);
    IHchange = [IHchange,temp];
end

nevokewidths = allwidths;
bad = isnan(allwidths) | ES==0 | isnan(IH);
allwidths = allwidths(~bad); ES = ES(~bad); IH = IH(~bad); dist = dist(~bad);
ESchange = ESchange(~bad); IHchange = IHchange(~bad);

lim = 12;

% Figure
s1 = subplot(5,2,3); histogram(nevokewidths,5:0.5:35,'normalization','probability','edgealpha',0)
hold on; histogram(allwidths,5:0.5:35,'normalization','probability','edgealpha',0)
xlabel('Spike width (a.u.)'); ylabel('Probability'); box off;
legend('All spikes','Evoked spikes'); legend boxoff;

bins = 0:0.01:0.4;
subplot(5,2,5); histogram(ES(allwidths>=lim),bins,'normalization','probability','edgealpha',0)
hold on; histogram(ES(allwidths<lim),bins,'normalization','probability','edgealpha',0)
xlim([0,0.4]); xlabel('Evoked spike probability'); ylabel('Probability');
legend('RS','FS'); legend boxoff; box off;
[h,p] = ttest2(ES(allwidths<lim),ES(allwidths>=lim))
[p,h] = ranksum(ES(allwidths<lim),ES(allwidths>=lim))

bins = 0:0.01:0.3;
s3 = subplot(5,2,6); histogram(IH(allwidths>=lim),bins,'normalization','probability','edgealpha',0)
hold on; histogram(IH(allwidths<lim),bins,'normalization','probability','edgealpha',0)
xlim([0,0.3]); xlabel('Inhibition intensity'); ylabel('Probability');
legend('RS','FS'); legend boxoff; box off;
[h,p] = ttest2(IH(allwidths<lim),IH(allwidths>=lim))
[p,h] = ranksum(IH(allwidths<lim),IH(allwidths>=lim))

bins = 0:0.5:4.5;
s2 = subplot(5,2,4); boxplot(dist,allwidths<lim,'notch','on','symbol','w.');
xticklabels({'RS','FS'}); ylabel('Distance (mm)');
set(findobj(gca,'type','line'),'linew',1.5);  box off;
xlim([0.5,2.5]); ylabel('Distance (mm)'); ylim([-0.2,4]);
[h,p] = ttest2(dist(allwidths<lim),dist(allwidths>=lim))
[p,h] = ranksum(dist(allwidths<lim),dist(allwidths>=lim))
nanmedian(dist(allwidths<lim))
nanmedian(dist(allwidths>=lim))
s2.Position(2:4) = s1.Position(2:4);
s2.Position(1) = s3.Position(1);

bins = 5:1:35;
subplot(5,2,7); histogram(allwidths(ESchange==0),bins,'normalization','probability','edgealpha',0)
hold on; histogram(allwidths(ESchange~=0),bins,'normalization','probability','edgealpha',0)
xlabel('Spike width (a.u.)'); ylabel('Probability'); 
legend('No change','ES change'); legend boxoff; box off;
subplot(5,2,8); histogram(allwidths(IHchange==0),bins,'normalization','probability','edgealpha',0)
hold on; histogram(allwidths(IHchange~=0),bins,'normalization','probability','edgealpha',0)
xlabel('Spike width (a.u.)'); ylabel('Probability'); 
legend('No change','IH change'); legend boxoff; box off;
subplot(5,2,9); histogram(allwidths(ESchange==-1),bins,'normalization','probability','edgealpha',0)
hold on; histogram(allwidths(ESchange==1),bins,'normalization','probability','edgealpha',0)
xlabel('Spike width (a.u.)'); ylabel('Probability'); 
legend('ES decrease','ES increase'); legend boxoff; box off;
subplot(5,2,10); histogram(allwidths(IHchange==-1),bins,'normalization','probability','edgealpha',0)
hold on; histogram(allwidths(IHchange==1),bins,'normalization','probability','edgealpha',0)
xlabel('Spike width (a.u.)'); ylabel('Probability'); 
legend('IH decrease','IH increase'); legend boxoff; box off;
set(findall(gcf,'-property','FontSize'),'FontSize',10)

%% Evoked vs Inhib
% Stimulation that evoked spikes vs stim that did not
evoked = []; 
nevoked = [];
for s = 1:length(SL)
    for c = 1:length(SL(s).chns)
        if(isempty(SL(s).ES))
            continue;
        end
        if(length(SL(s).ES{c})/length(SL(s).stim) < 0.05)
            continue;
        end
        temp = 1:length(SL(s).stim);
        inds = discretize(SL(s).ES{c},SL(s).stim);
        inds(inds>max(temp)) = []; inds(isnan(inds)) = [];
        temp(inds) = [];
        evoked(end+1) = nanmedian(SL(s).IH{c}(inds));
        nevoked(end+1) = nanmedian(SL(s).IH{c}(temp));
    end
end
bad = isnan(evoked) | isnan(nevoked);
evoked = evoked(~bad); nevoked = nevoked(~bad);
figure; boxplot([evoked,nevoked]*1000,[ones(1,length(evoked)),2*ones(1,length(nevoked))],...
    'notch','on','symbol','w');
ylim([0,130])
set(findobj(gca,'type','line'),'linew',2); 
xticklabels({'Evoked','Not evoked'});
ylabel('Inhibition (ms)');
[h,p] = ttest(evoked,nevoked)

% All combined
ESslope = [];
IHslope = [];
Dist = [];
ESper = [];
EStime = [];
ESpval = [];
IHpval = [];
minp = 0.01;
mindur = 600;
for s = 1:length(SL)
    
    if(SL(s).stim(end)-SL(s).stim(1) < mindur)
        continue;
    end
    
    EStemp = SL(s).ESslope; EStemp(SL(s).ESpval>=minp) = 0;
    ESslope = [ESslope,sign(EStemp)];
    IHtemp = SL(s).IHslope; IHtemp(SL(s).IHpval>=minp) = 0;
    IHslope = [IHslope,sign(IHtemp)];
    Dist = [Dist,SL(s).dist];
    
    temp = cellfun(@(x) length(x)/length(SL(s).stim),SL(s).ES);
    ESper = [ESper,temp];
    
    temp = cellfun(@nanmedian, SL(s).ES_delay);
    EStime = [EStime,temp];
    
end

bad = isnan(ESslope) | isnan(IHslope);% | ESper<0.03;
bad = bad | ESper<0.2;
ESslope(bad) = [];
IHslope(bad) = [];
Dist(bad) = [];
ESper(bad) = [];
EStime(bad) = [];

% By distance. Good! 
mean(Dist((ESslope==-1 & IHslope ==-1) | (ESslope==1 & IHslope ==1)))
mean(Dist((ESslope==1 & IHslope ==-1) | (ESslope==-1 & IHslope ==1)))

groups = ones(1,length(Dist));
groups((ESslope==1 & IHslope ==-1) | (ESslope==-1 & IHslope ==1) ) = 2;
groups((ESslope==0 | IHslope ==0)) = 3;

figure; 
subplot(1,3,1); boxplot(Dist,groups,'notch','on','symbol','w');
xticklabels({'Covary','Antivary','Neither'})
ylabel('Distance (mm)');
set(findobj(gca,'type','line'),'linew',2); 
xl = xlim; avg = median(Dist); hold on; plot(xl,[avg,avg],'k--');
p = [];
[~,p(1)] = ttest2(Dist(groups==1),Dist(groups==2));
[~,p(2)] = ttest2(Dist(groups==2),Dist(groups==3));
[~,p(3)] = ttest2(Dist(groups==1),Dist(groups==3));
sigstar({[1,2],[2,3],[1,3]},p,0,20,0);

% ES percentage
subplot(1,3,3) ; boxplot(ESper*100,groups,'notch','on','symbol','w');
xticklabels({'Covary','Antivary','Neither'})
ylabel('Probability (%)');
set(findobj(gca,'type','line'),'linew',2); 
xl = xlim; avg = nanmedian(ESper*100); hold on; plot(xl,[avg,avg],'k--');
ylim([-2,70])
p = [];
[~,p(1)] = ttest2(ESper(groups==1),ESper(groups==2));
[~,p(2)] = ttest2(ESper(groups==2),ESper(groups==3));
[~,p(3)] = ttest2(ESper(groups==1),ESper(groups==3));
sigstar({[1,2],[2,3],[1,3]},p,0,20,0);

% ES timing
subplot(1,3,2); boxplot(EStime*1000,groups,'notch','on','symbol','w');
xticklabels({'Covary','Antivary','Neither'})
ylabel('Time (ms)');
set(findobj(gca,'type','line'),'linew',2); 
xl = xlim; avg = median(EStime*1000); hold on; plot(xl,[avg,avg],'k--');
ylim([0,20])
p = [];
[~,p(1)] = ttest2(EStime(groups==1),EStime(groups==2));
[~,p(2)] = ttest2(EStime(groups==2),EStime(groups==3));
[~,p(3)] = ttest2(EStime(groups==1),EStime(groups==3));
sigstar({[1,2],[2,3],[1,3]},p,0,20,0);

% For checking distribution
figure; 
subplot(1,3,1);
histogram(Dist(groups==1),0:0.4:4,'normalization','probability'); 
hold on; histogram(Dist(groups==2),0:0.4:4,'normalization','probability'); 
% hold on; histogram(Dist(groups==3),0:0.4:4,'normalization','probability'); 
xlabel('Distance (mm)'); ylabel('Proportion of spikes');
legend('Covary','Antivary'); legend('box','off');

subplot(1,3,3);
histogram(ESper(groups==1)*100,0:5:100,'normalization','probability'); 
hold on; histogram(ESper(groups==2)*100,0:5:100,'normalization','probability'); 
xlabel('Probability (%)'); ylabel('Proportion of spikes');

subplot(1,3,2);
histogram(EStime(groups==1),0:0.001:0.02,'normalization','probability'); 
hold on; histogram(EStime(groups==2),0:0.001:0.02,'normalization','probability'); 
xlabel('Time (ms)'); ylabel('Proportion of spikes');

%% Magnitude of changes
ESslope = [];
IHslope = [];
ES = [];
IH = [];
dur = [];
for s = 1:length(SL)
    ESslope = [ESslope,SL(s).ESslope];
    IHslope = [IHslope,SL(s).IHslope];
    
    ES = [ES, cellfun(@length, SL(s).ES)./length(SL(s).stim)];
    IH = [IH, cellfun(@nanmedian, SL(s).IH)];
    
    dur = [dur,SL(s).stim(end)-SL(s).stim(1)];
end

bad = ES < 0.05;
ESslope(bad) = [];
IHslope(bad) = [];
ES(bad) = [];
IH(bad) = [];


chnind = {};
mindur = 600;
esper = 0.2;
minpval = 0.01;
for s = 1:length(SL)
    if(SL(s).stim(end)-SL(s).stim(1)<mindur)
        continue;
    end
    
    temp1 = cellfun(@length, SL(s).ES)./length(SL(s).stim);
    temp2 = SL(s).ESpval;
    
    chnind{s} = find(temp1>=esper & temp2<=minpval);
    
end

%% Covarying
% Interesting. Yes relationship for evoked spike changes, no for
% inhibition. More proof the two are independent? 
cov = [];
ncov = [];
icov = [];
incov = [];
allcov = 0;
allcovmatch = 0;
allpairs = 0;
for s = 1:length(SL)
    covmatch = 0;
    icovmatch = 0;
    ncovmatch = 0;
    incovmatch = 0;
    covtot = 0;
    ncovtot = 0;
    for i = 1:(length(SL(s).chns)-1)
        if((length(SL(s).ES{i})/length(SL(s).stim))<0.05)
            continue;
        end
        for j = (i+1):length(SL(s).chns)
            if((length(SL(s).ES{j})/length(SL(s).stim))<0.05)
                continue;
            end
            allpairs = allpairs+1;
            if(SL(s).Cov(i,j)==1)
                allcov = allcov+1;
                covtot = covtot+1;
                if sign(SL(s).ESslope(i)) == sign(SL(s).ESslope(j))
                    covmatch = covmatch+1;
                    allcovmatch = allcovmatch+1;
                end
                if sign(SL(s).IHslope(i)) == sign(SL(s).IHslope(j))
                    icovmatch = icovmatch+1;
                end
            elseif (SL(s).Cov(i,j)==-1)
                ncovtot = ncovtot+1;
                if sign(SL(s).ESslope(i)) == sign(SL(s).ESslope(j))
                    ncovmatch = ncovmatch+1;
                end
                if sign(SL(s).IHslope(i)) == sign(SL(s).IHslope(j))
                    incovmatch = incovmatch+1;
                end
            end
        end
    end
    cov = [cov,covmatch/covtot];
    ncov = [ncov,ncovmatch/ncovtot];
    icov = [icov,icovmatch/covtot];
    incov = [incov,incovmatch/ncovtot];
end

figure; 
subplot(1,2,1);
boxplot([cov,ncov],[ones(1,length(cov)),2*ones(1,length(ncov))],'notch','on',...
    'symbol','w');
ylabel('Probabillity of same change');
xticklabels({'Covary','Not Covary'}); ylim([0,1])
set(findobj(gca,'type','line'),'linew',2); 
title('Evoked Spike');
[p,~] = ranksum(cov,ncov);
if(p<=0.05)
    sigstar([1,2],p,0,20,0);
end

subplot(1,2,2);
boxplot([icov,incov],[ones(1,length(icov)),2*ones(1,length(incov))],'notch','on',...
    'symbol','w');
ylabel('Probabillity of same change');
xticklabels({'Covary','Not Covary'}); ylim([0,1])
set(findobj(gca,'type','line'),'linew',2); 
title('Inhibition');
[p,~] = ranksum(icov,incov);
sigstar([1,2],p,0,20,0);

ESslope = [];
IHslope = [];

for s = 1:length(SL)
    ESslope = [ESslope,SL(s).ESslope];
    IHslope = [IHslope,SL(s).IHslope];

end

%% Distance and same change
bins = 0.25:0.25:3;
binned = {};
ibinned = {};
for s = 1:length(SL)
    chns = SL(s).chns;
    pos = [];
    for i = 1:length(chns)
        [pos(i,1),pos(i,2),~] = GetWadeChannelPosition(chns(i));
    end
    
    dist = pdist2(pos,pos)*0.4;
    dist(dist==0) = -1;
    
    binind = discretize(dist,bins);
    data = cell(length(bins)-1,1);
    
    for i = 1:(length(chns)-1)
        if isnan(SL(s).ESslope(i)) || sign(SL(s).ESslope(i))==0
            continue;
        end
        for j = i:length(chns)
            if isnan(SL(s).ESslope(j))
                continue;
            end
            if(isnan(binind(i,j)))
                continue;
            end
            same = sign(SL(s).ESslope(i)) == sign(SL(s).ESslope(j));
            data{binind(i,j)} = [data{binind(i,j)},same];
            
        end
    end
    
    binned{s} = data;
    
    idata = cell(length(bins)-1,1);
    
    for i = 1:(length(chns)-1)
        if isnan(SL(s).IHslope(i)) || sign(SL(s).IHslope(i))==0
            continue;
        end
        for j = i:length(chns)
            if isnan(SL(s).IHslope(j))
                continue;
            end
            
            if(isnan(binind(i,j)))
                continue;
            end
            
            same = sign(SL(s).IHslope(i)) == sign(SL(s).IHslope(j));
            idata{binind(i,j)} = [idata{binind(i,j)},same];
            
        end
    end
    
    ibinned{s} = idata;
end

ratio = cellfun(@(x) cellfun(@(y) sum(y==1)./length(y),x),binned,'uniformoutput',false);
ratio = cell2mat(ratio);

iratio = cellfun(@(x) cellfun(@(y) sum(y==1)./length(y),x),ibinned,'uniformoutput',false);
iratio = cell2mat(iratio);

figure; 
subplot(1,2,1); boxplot(ratio','notch','on','symbol','w');
xticks(1.5:2:(length(bins)));
xticklabels(bins(2:2:end)); xtickangle(45);
xlabel('Distance (mm)'); ylabel('Probability of same change'); 
title('Evoked Spikes'); set(findobj(gca,'type','line'),'linew',2); 
xlim([0.5,11.5]); ylim([-0.05,1.05]);


subplot(1,2,2); boxplot(iratio','notch','on','symbol','w');
xticks(1.5:2:(length(bins)));
xticklabels(bins(2:2:end));xtickangle(45);
xlabel('Distance (mm)'); ylabel('Probability of same change'); 
title('Inhibition'); set(findobj(gca,'type','line'),'linew',2); 
xlim([0.5,11.5]); ylim([-0.05,1.05]);

%% Espikes vs most recent spike 
bw = 0.001;
coefs = [];
total = 0;
allrates = [];
spikenum = [];
dur = [];
pval = [];
track = [];
for s = 1:length(SL)
    if(SL(s).stimfreq > 10)
        continue;
    end
    
    stim = SL(s).stim;
    dur = [dur,max(stim)-min(stim)];
    if(dur(end)<1200)
        continue;
    end
    maxbin = min(0.05,1/SL(s).stimfreq);
    bins = 0.005:bw:maxbin;
    
    for i = 1:length(SL(s).ES)
        if(isempty(SL(s).ES{i}))
            continue;
        end
        
        ESper = length(SL(s).ES{i})/length(stim);
        if(ESper < 0.05 || ESper > 0.7)
            continue;
        end
        
        spikes = SL(s).spiketimes{i};
        if(isempty(spikes)), continue; end
        spikerate = length(spikes)/(spikes(end)-spikes(1));
%         if(spikerate<10 || spikerate > 30)
%             continue;
%         end
        
        ES = SL(s).ES{i};
        allspikes = sort([0;ES;spikes;inf]);
        
        ESstim = discretize(ES,stim);
        ESstim = unique(ESstim); ESstim(isnan(ESstim)) = [];
        
        allrecent = discretize(stim,allspikes);
        allrecent = stim'-allspikes(allrecent);
        
        ESrecent = discretize(stim(ESstim),allspikes);
        ESrecent = stim(ESstim)'-allspikes(ESrecent);
        
        allcounts = histcounts(allrecent,bins);
        EScounts = histcounts(ESrecent,bins);
        ratio = EScounts./allcounts;
        
        if(any(allcounts < 100) || any(ratio<0.15))
            continue;
        end
        
        if(sum(EScounts<20)>10)
            continue;
        end
        
        ratio = smooth(ratio,5,'rlowess');
        
        [cor, lags] = CrossCorr(spikes, 'ts2', spikes, 'binsize', bw, 'lag', [min(bins) max(bins)]-bw/2);
        
        cor = smooth(cor,5,'rlowess');
        
        bad = isnan(ratio);
        if(sum(bad)>(length(ratio)/2))
            continue;
        end

        
        [R,p] = corrcoef(ratio(~bad),cor(~bad));
        total = total+1;
        pval = [pval;p(1,2)];
        track(end+1,:) = [s,i];
        coefs = [coefs;R(1,2)];
        allrates = [allrates,spikerate];

    end
    
end

figure; histogram(coefs);
figure; boxplot(coefs,'notch','on','symbol','w');

%% Spikes per stim

s = SL(27);
ES = s.ES; ES = ES(cellfun(@(x) ~isempty(x),ES));
ESmat = zeros(length(ES),length(s.stim));
PrevSpike = zeros(length(ES),length(s.stim));
for e = 1:length(ES)
    ind = discretize(ES{e},s.stim);
    ind(isnan(ind)) = [];
    ESmat(e,unique(ind)) = 1;
    
    bins = [floor(s.stim(1));ES{e};inf]';
    temp = discretize(s.stim,bins);
    
    PrevSpike(e,:) = s.stim-bins(temp);
end

keep = sum(ESmat,2)/length(s.stim) > 0.1 & sum(ESmat,2)/length(s.stim) < 0.5;

PrevSpike = log(PrevSpike(keep,:));
PrevSpike = zscore(PrevSpike,[],2);
ESmat = ESmat(keep,:);


temp = sum(ESmat);
[~,ind] = sort(temp);


[coeff,score,latent,tsquared,explained,mu] = pca(PrevSpike');%,'Algorithm','eig');
variance = cumsum(explained/sum(explained));

figure; scatter3(score(:,1),score(:,2),score(:,3),'.');

temp = logical(ESmat(6,:));
figure; scatter3(score(temp,1),score(temp,2),score(temp,3),'.');
hold on; scatter3(score(~temp,1),score(~temp,2),score(~temp,3),'.');

%% Evoked vs inhib with correlation
Dist = [];
ESper = [];
EStime = [];
cor = [];
bw = 1; 
for s = 1:length(SL)
    bins = SL(s).stim(1):bw:SL(s).stim(end);  
    stimcount = histcounts(SL(s).stim,bins);
    stimind = discretize(SL(s).stim,bins);
    disp(s)
    for c  = 1:length(SL(s).ES)
        if isempty(SL(s).ES{c})
            continue;
        end
        
        evokedcount = histcounts(SL(s).ES{c},bins);
        ES = evokedcount./stimcount;
        ES(isinf(ES)) = nan;
        
        ES = smooth(ES,30);
        
        IH = SL(s).IH{c};
        IH = arrayfun(@(x) nanmean(IH(stimind==x)),1:(length(bins)-1));
        
        IH = smooth(IH,30);
        
        [R,P] = corrcoef(ES,IH);
        if(P(1,2)<0.05)
            cor(end+1) = sign(R(1,2));
        else
            cor(end+1) = 0;
        end
        
        Dist(end+1) = SL(s).dist(c);
        ESper(end+1) = length(SL(s).ES{c})/length(SL(s).stim);
        EStime(end+1) = mean(SL(s).ES_delay{c});
        
    end
    
end

bad = ESper<0.03;

groups = zeros(1,length(Dist));
groups(cor==1) = 2;
groups(cor==-1) = 3;
groups(cor==0) = 1;

% Distance
figure; 
subplot(1,2,1); boxplot(Dist(~bad),groups(~bad),'notch','on','symbol','w');
xticklabels({'Uncorr','Corr','Anti'})
ylabel('Distance (mm)');
set(findobj(gca,'type','line'),'linew',2); 
% xl = xlim; avg = median(Dist); hold on; plot(xl,[avg,avg],'k--');
p = [];
[p(1),~] = ranksum(Dist(groups==1),Dist(groups==2));
[p(2),~] = ranksum(Dist(groups==2),Dist(groups==3));
[p(3),~] = ranksum(Dist(groups==1),Dist(groups==3));
sigstar({[2,3]},p(2),0,20,0);
box off;


% ES percentage
subplot(1,3,2) ; boxplot(ESper(~bad)*100,groups(~bad),'notch','on','symbol','w');
xticklabels({'Corr','Anti','Neither'})
ylabel('Probability (%)');
set(findobj(gca,'type','line'),'linew',2); 
% xl = xlim; avg = nanmedian(ESper(~bad)*100); hold on; plot(xl,[avg,avg],'k--');
ylim([-2,60])
p = [];
[p(1),~] = ranksum(ESper(~bad & groups==1),ESper(~bad & groups==2));
[p(2),~] = ranksum(ESper(~bad & groups==2),ESper(~bad & groups==3));
[p(3),~] = ranksum(ESper(~bad & groups==1),ESper(~bad & groups==3));
sigstar({[1,2],[1,3]},p([1,3]),0,20,0);
box off;


% ES timing
subplot(1,3,3); boxplot(EStime(~bad)*1000,groups(~bad),'notch','on','symbol','w');
xticklabels({'Corr','Anti','Neither'})
ylabel('Evoked Spike Time (ms)');
set(findobj(gca,'type','line'),'linew',2); 
% xl = xlim; avg = median(EStime*1000); hold on; plot(xl,[avg,avg],'k--');
ylim([0,18])
p = [];
[p(1),~] = ranksum(EStime(groups==1),EStime(groups==2));
[p(2),~] = ranksum(EStime(groups==2),EStime(groups==3));
[p(3),~] = ranksum(EStime(groups==1),EStime(groups==3));
% sigstar({[1,2],[2,3],[1,3]},p,0,20,0);
box off;


subplot(2,3,4);
histogram(Dist(groups==1&~bad),0:0.4:4,'normalization','probability'); 
hold on; histogram(Dist(groups==2&~bad),0:0.4:4,'normalization','probability'); 
% hold on; histogram(Dist(groups==3),0:0.4:4,'normalization','probability'); 
xlabel('Distance (mm)'); ylabel('Proportion of spikes');
legend('Correlated','Anticorrelated'); legend('box','off');
box off;

subplot(2,3,5);
histogram(EStime(groups==1&~bad),0:0.001:0.02,'normalization','probability'); 
hold on; histogram(EStime(groups==2&~bad),0:0.001:0.02,'normalization','probability'); 
xlabel('Time (ms)'); ylabel('Proportion of spikes');
box off;

subplot(2,3,6);
histogram(ESper(groups==1)*100,0:5:100,'normalization','probability'); 
hold on; histogram(ESper(groups==2)*100,0:5:100,'normalization','probability'); 
xlabel('Probability (%)'); ylabel('Proportion of spikes');
box off;

% Normalized
bins = 0:0.4:4;
total = histcounts(Dist(~bad),'BinEdges',bins);
PosCor = histcounts(Dist(groups==1&~bad),'BinEdges',bins);
NegCor = histcounts(Dist(groups==2&~bad),'BinEdges',bins);
Neither = histcounts(Dist(groups==3&~bad),'BinEdges',bins);

figure; plot(bins(1:end-1)+0.2,PosCor./total); 
hold on; plot(bins(1:end-1)+0.2,NegCor./total);
hold on; plot(bins(1:end-1)+0.2,Neither./total);
legend('PosCor','NegCor','Neither');
xlabel('mm'); ylabel('Proportion of all spikes');

%% ES vs IH with FR & ES correlated spikes removed
Dist = [];
ESper = [];
EStime = [];
cor = [];
rval = [];
bw = 1; 
for s = 1:length(SL)
    bins = SL(s).stim(1):bw:SL(s).stim(end);  
    stimcount = histcounts(SL(s).stim,bins);
    stimind = discretize(SL(s).stim,bins);
    disp(s)   
    for c  = 1:length(SL(s).ES)
        if isempty(SL(s).ES{c})
            continue;
        end
        
        evokedcount = histcounts(SL(s).ES{c},bins);
        ES = evokedcount./stimcount;
        ES(isinf(ES)) = nan;
        
        ES = smooth(ES,30);
        
        FR = histcounts(SL(s).spiketimes{c},bins);
        [R,P] = corrcoef(ES,FR);
        
        if(P(1,2)<0.05), continue; end
        
        IH = SL(s).IH{c};
        IH = arrayfun(@(x) nanmean(IH(stimind==x)),1:(length(bins)-1));
        
        IH = smooth(IH,30);
        
        [R,P] = corrcoef(IH,FR);
        
        if(P(1,2)<0.05), continue; end
        
        [R,P] = corrcoef(ES,IH);
        if(P(1,2)<0.05)
            cor(end+1) = sign(R(1,2));
        else
            cor(end+1) = 0;
        end
        rval(end+1) = R(1,2);
        
        Dist(end+1) = SL(s).dist(c);
        ESper(end+1) = length(SL(s).ES{c})/length(SL(s).stim);
        EStime(end+1) = mean(SL(s).ES_delay{c});
        
    end
    
end

bad = ESper<0.03;

groups = zeros(1,length(Dist));
groups(cor==1) = 2;
groups(cor==-1) = 3;
groups(cor==0) = 1;

% Distance
subplot(1,2,2);
boxplot(Dist(~bad),groups(~bad),'notch','on','symbol','w');
xticklabels({'Uncorr','Corr','Anti'})
ylabel('Distance (mm)');
set(findobj(gca,'type','line'),'linew',2); 
% xl = xlim; avg = median(Dist); hold on; plot(xl,[avg,avg],'k--');
p = [];
[p(1),~] = ranksum(Dist(groups==1),Dist(groups==2));
[p(2),~] = ranksum(Dist(groups==2),Dist(groups==3));
[p(3),~] = ranksum(Dist(groups==1),Dist(groups==3));
box off;

% R values
bad = cor==0;
tempr = rval; tempr(bad) = []; tempr = abs(tempr);
tempd = Dist; tempd(bad) = [];

bins = 0.25:0.5:3;
ind = discretize(tempd,bins);
avg = arrayfun(@(x) nanmedian(tempr(ind==x)), 1:max(ind));

figure; 
plot(avg,'color',[0.4,0.4,0.4]); hold on;
boxplot(tempr,ind,'notch','on','symbol','w'); 
set(findobj(gca,'type','line'),'linew',1.5); 
ylabel('R value'); xlabel('Distance (mm)')
xticks(1.5:2:length(bins));
xticklabels(bins(2:2:length(bins))); 
box off;

%% Changes wrt distance from stimulated site
set(0,'defaultAxesFontSize',12)

D = [];
ES = [];
IH = [];
problim = 0.03;
for s = 1:length(SL)
    ind = find(cellfun(@(x) ~isempty(x), SL(s).ES));
    temp = SL(s).ES(ind);

    prob = cellfun(@(x) length(x)/length(SL(s).stim),temp);
    ind(prob<problim) = [];
    
    D = [D,SL(s).dist(ind)];
    es = sign(SL(s).ESslope(ind));
    es(SL(s).ESpval > 0.05) = 0;
    ES = [ES,es];
    
    ih = sign(SL(s).IHslope(ind));
    ih(SL(s).IHpval > 0.05) = 0;
    IH = [IH,ih];
    
end

bins = 0.25:0.5:3;
ind = discretize(D,bins);

ESchange = arrayfun(@(x) sum(ES(ind==x)~=0)/sum(ind==x), 1:max(ind));
ESdecrease = arrayfun(@(x) sum(ES(ind==x)==-1)/sum(ES(ind==x)~=0), 1:max(ind));
subplot(2,2,3); plot(bins(1:end-1)+mean(diff(bins))/2,ESchange*100); 
hold on; plot(bins(1:end-1)+mean(diff(bins))/2,ESdecrease*100);
xlabel('Distance (mm)'); ylabel('Probability (%)'); box off;
xlim([0.25,2.75])

IHchange = arrayfun(@(x) sum(IH(ind==x)~=0)/sum(ind==x), 1:max(ind));
IHdecrease = arrayfun(@(x) sum(IH(ind==x)==-1)/sum(IH(ind==x)~=0), 1:max(ind));
subplot(2,2,4); plot(bins(1:end-1)+mean(diff(bins))/2,IHchange*100); 
hold on; plot(bins(1:end-1)+mean(diff(bins))/2,IHdecrease*100);
xlabel('Distance (mm)'); box off; xlim([0.25,2.75])

% Sums
val = arrayfun(@(x) sum(IH(ind==x)==1), 1:max(ind)); val = val(:);


% Stats
ESchangepvals = [];
ESdecreasepvals = [];
IHchangepvals = [];
IHdecreasepvals = [];
for i = 1:max(ind)
    for j = 1:max(ind)
        n1 = sum(ES(ind==i) ~= 0); n2 = sum(ES(ind==i) == 0);
        n3 = sum(ES(ind==j) ~= 0); n4 = sum(ES(ind==j) == 0);
        
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,ESchangepvals(i,j)] = crosstab(x1,x2);
        
        n1 = sum(ES(ind==i) == -1); n2 = sum(ES(ind==i) == 1);
        n3 = sum(ES(ind==j) == -1); n4 = sum(ES(ind==j) == 1);
        
        [h,decreasepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,ESdecreasepvals(i,j)] = crosstab(x1,x2);
        
        
        n1 = sum(IH(ind==i) ~= 0); n2 = sum(IH(ind==i) == 0);
        n3 = sum(IH(ind==j) ~= 0); n4 = sum(IH(ind==j) == 0);
        
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,IHchangepvals(i,j)] = crosstab(x1,x2);
        
        n1 = sum(IH(ind==i) == -1); n2 = sum(IH(ind==i) == 1);
        n3 = sum(IH(ind==j) == -1); n4 = sum(IH(ind==j) == 1);
        
        [h,decreasepvalsf(i,j)] = fishertest([n1,n2; n3,n4]);
        x1 = [repmat('a',n1+n2,1); repmat('b',n3+n4,1)];
        x2 = [repmat(1,n1,1); repmat(2,n2,1); repmat(1,n3,1); repmat(2,n4,1)];
        [tbl,chi,IHdecreasepvals(i,j)] = crosstab(x1,x2);
        
    end
end

%% Firing rate
ESper = [];
ES = {};
FR = {};
IH = {};
Cor = [];
ICor = [];
bw = 1; 
ESP = [];
IHP = [];
D = [];
for s = 1:length(SL)
    bins = SL(s).stim(1):bw:SL(s).stim(end);  
    stimcount = histcounts(SL(s).stim,bins);
    stimind = discretize(SL(s).stim,bins);
    disp(s)
    for c  = 1:length(SL(s).ES)
        if isempty(SL(s).ES{c})
            continue;
        end
        
        evokedcount = histcounts(SL(s).ES{c},bins);
        temp = evokedcount./stimcount;
        temp(isinf(temp)) = 0;
        ES{end+1} = temp;
        
        FR{end+1} = histcounts(SL(s).spiketimes{c},bins);
        
        ESper(end+1) = length(SL(s).ES{c})/length(SL(s).stim);
               
        [R,P] = corrcoef(ES{end},FR{end});
               
        if(ESper(end) < 0.03)
            Cor(end+1) = nan;
        elseif(P(1,2)<0.05)
            Cor(end+1) = sign(R(1,2));
        else
            Cor(end+1) = 0;
        end
        ESP(end+1) = P(1,2);
        
        IH{end+1} = arrayfun(@(x) nanmedian(SL(s).IH{c}(stimind==x)),1:max(stimind));

        [R,P] = corrcoef(IH{end},FR{end}(1:length(IH{end})));
               
        if(ESper(end) < 0.03)
            ICor(end+1) = nan;
        elseif(P(1,2)<0.05)
            ICor(end+1) = sign(R(1,2));
        else
            ICor(end+1) = 0;
        end
        IHP(end+1) = P(1,2);
        
        D(end+1) = SL(s).dist(c);
        
    end
    
end

% Figures
avgFR = cellfun(@nanmean, FR);
sig = find(ESP<0.001 & IHP<0.001);
[~,sorted] = sort(avgFR(sig));
ind = sig(sorted(end));

figure;
subplot(3,1,1); colororder({'k','b'})
plot(smooth(FR{ind},20),'linewidth',1.5); 
ylabel('Firing rate (Hz)');
yyaxis right; plot(smooth(ES{ind},20)*100,'linewidth',1.5); 
xticklabels({}); xlim([0,600]); set(gca,'FontSize',11); box off;
ylabel('Evoked Spike Probability (%)');
subplot(3,1,2); set(gca,'colororder',[0 0 0; 0 0 1])
plot(smooth(FR{ind},20),'linewidth',1.5); 
ylabel('Firing rate (Hz)');
yyaxis right; plot(smooth(IH{ind},20)*1000,'linewidth',1.5); 
xlim([0,600]); set(gca,'FontSize',11); 
xlabel('Time (s)'); ylabel('Inhibition Duration (ms)'); box off;

% Correlated vs anticorrelated
bad = isnan(Cor);
dtemp = D; dtemp(bad) = [];
temp = Cor; temp(bad) = [];
ind = zeros(1,length(temp)); ind(temp==1) = 1; ind(temp==-1) = 2;
figure('Position',[680, 300, 1000, 500]); 
subplot(1,2,1); boxplot(dtemp,ind,'notch','on','symbol','w');
[p(1),~] = ranksum(dtemp(temp==0),dtemp(temp==1));
[p(2),~] = ranksum(dtemp(temp==0),dtemp(temp==-1));
[p(3),~] = ranksum(dtemp(temp==1),dtemp(temp==-1));
set(findobj(gca,'type','line'),'linew',2);
sigstar({[1,2],[1,3],[2,3]},p,0,20,0);
xticklabels({'Uncorrelated','Correlated','Anticorrelated'})
ylabel('Distance (mm)'); yl = ylim;
title('Evoked Spikes and Firing Rate');

% IH
subplot(1,2,2);
bad = isnan(IHP);
dtemp = D; dtemp(bad) = [];
temp = IHP; temp(bad) = [];
boxplot(dtemp,temp>=0.05,'notch','on','symbol','w');
[p,~] = ranksum(dtemp(temp<0.05),dtemp(temp>=0.05));
set(findobj(gca,'type','line'),'linew',2);
xticklabels({'Uncorrelated','Anticorrelated'}); ylim(yl);
title('Inhibition and Firing Rate')

%% Latency change over time and stim rate
Dist = [];
Freq = [];
PVals = [];
Slopes = [];
perlim = 0.03;
for s = 1:length(SL)
    disp(s);
    freq = length(SL(s).stim)/(SL(s).stim(end)-SL(s).stim(1));
    for c = 1:length(SL(s).ES_delay)
        if(length(SL(s).ES{c})<2), continue; end
        ESper = length(SL(s).ES{c})/length(SL(s).stim);
        if(ESper < perlim || ESper > 0.95), continue; end
        
        [slope,pval] = extractTrend(SL(s).ES_delay(c),SL(s).ES{c},0);
        PVals(end+1) = pval;
        Slopes(end+1) = slope;
        Freq(end+1) = freq;
        Dist(end+1) = SL(s).dist(c);
    end
end

figure; subplot(1,2,1);
boxplot(Dist,isnan(PVals), 'notch','on','symbol','w');
[p,~] = ranksum(Dist(isnan(PVals)),Dist(~isnan(PVals)));
sigstar({[1,2]},p,0,20,0);
xticklabels({'Change','No Change'}); title('Latency changes and distance');
ylabel('Distance (mm)');
set(findobj(gca,'type','line'),'linew',1.5);

subplot(1,2,2);
boxplot(Freq,isnan(PVals), 'notch','on');
[p,~] = ranksum(Freq(isnan(PVals)),Freq(~isnan(PVals)));
sigstar({[1,2]},p,0,20,0);
xticklabels({'Change','No Change'}); title('Latency changes and stim frequency');
ylabel('Stim frequency (Hz)');
set(findobj(gca,'type','line'),'linew',1.5);


pos = find(sign(Slopes) == 1 & PVals<0.05);
neg = find(sign(Slopes) == -1 & PVals<0.05);







