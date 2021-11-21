%% Sorting all spikes offline to ensure online sorting matches

%% Parse all data
clear; close all;

%% Load
% Path all folders are in
basepath = 'R:\Yun\Kronk\Ripple';

days = dir(basepath);

% Remove parent directory
days = days(3:end);

% Needs to be a subdirectory
isdir = cell2mat(extractfield(days,'isdir'));

days = days(isdir);

% Needs to be numeric
name = extractfield(days,'name');
days = name(cellfun(@(x) ~isnan(str2double(x)),name));

verbose = 1; %1 to print all times, 0 to not.

for d = 1:length(days)
    
    path = fullfile(basepath,days{d});
    files = dir(fullfile(path,'*.nev')); %find all ripple files
    files = extractfield(files,'name');
    
    for f = 1:length(files)
        
        try % just to get the ones that'll work
            
            %% Check to see if relevant file
            if(~(contains(files{f},'5Hz') || contains(files{f},'10Hz') || contains(files{f},'20Hz') ||...
                    contains(files{f},'Gamm')))
                continue;
            end
            
            filename = files{f}(1:end-4);
            
            if(exist(fullfile(path,[filename,'_OfflineSort.mat'])))
                continue;
            end
            
            %% Check to see if SD exists and load
            disp([num2str(d),'. ',days{d},' - ',num2str(f),'. ',files{f}]);
            
            if(verbose), tic; fprintf('Loading SD and spikes...'); end
            
            sdfile = fullfile(path,[filename,'_SD.mat']);
            if(~exist(sdfile))
                continue;
            else
                load(sdfile);
            end
            
            %% Load spikes (have to for the ones that are offset)
            spikefile = fullfile(path,[filename,'_Spikes.mat']);
            if(exist(spikefile))
                load(spikefile)
            else
                [~, spikes, ~, ~, ~, ~, ~, ~, ~, ~,~]...
                    = LoadRippleData('path',path,'filename',filename,'flag',2);
                save(spikefile,'spikes');
            end
            [chns,codes,~] = findSpikeChnCode(spikes);
            
            if(verbose), t = toc; fprintf(' %f seconds \n', t); end
            
            %% Loop through channels
            % Set some variables
            fs = 30000;
            
            spikerange = round(-0.0005*fs):1:round(0.0015*fs);
            
            sortedspikes = {};
            spikeshape = {};
            allchns = unique(chns);
            for c = 1:length(allchns)
                
                %% Load initial data and clean
                fprintf(' Channel %d \n',allchns(c));
                
                
                nspikes = sum(chns==allchns(c));
                
                % Get threshold from previous traces
                chnind = find(chns==allchns(c));
                
                if(all(cellfun(@isempty,SD.traces(chnind))))
                    sortedspikes(chnind) = cell(1,length(chnind));
                    spikeshape(chnind) = cell(1,length(chnind));
                    continue;
                end
                
                trough = cellfun(@(x) min(x(5:30)), SD.traces(chnind),'uniformoutput',false);
                
                bad = cellfun(@(x) isempty(x)|| isnan(x) || x<-500,trough);
                % Skip if no good spike traces
                if(sum(bad) == length(trough))
                    sortedspikes(chnind) = cell(1,length(chnind));
                    spikeshape(chnind) = cell(1,length(chnind));
                    continue;
                end
                trough = cell2mat(trough(~bad));
                threshold = 0.75*min(trough);
                
                %% Load first 5 minutes to do initial sorting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(verbose), tic; fprintf('  Loading initial data...'); end
                
                [wave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~, session_time] = ...
                    LoadRippleData('path',path,'filename',filename,'flag',1,...
                    'wavechannel',allchns(c),'times',[0,300]);
                
                % Filter data
                temp = HPF(wave,fs,916);
                
                [snips, cross] = GetSnips(temp, threshold, spikerange, 1000, SD.stim,  0.0012, fs);
                
                if(size(snips,2) < 10)
                    sortedspikes(chnind) = cell(1,length(chnind));
                    spikeshape(chnind) = cell(1,length(chnind));
                    continue;
                end
                
                if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                
                %% 2 window discrim if only 1 spike
                % For 2 window discrim
                spiketimes = cell(1,length(chnind));
                if(nspikes==1)
                    
                    if(verbose), tic; fprintf('  2 window discrim...'); end
                    
                    % Use snips that don't have any nan and are >2ms from stim time for parameters
                    ind = discretize(cross/fs,SD.stim);
                    bad = find(~isnan(ind));
                    lag = cross(bad)/fs-SD.stim(ind(bad));
                    bad = bad(lag<0.002);
                    
                    goodsnips = snips; goodsnips(:,bad) = [];
                    goodsnips(:,any(isnan(goodsnips))) = [];
                    
                    avg = nanmean(goodsnips,2);
                    
                    % With 2 window discrim
                    % Get parameters
                    p = struct;
                    [~,minind] = min(avg);
                    [~,maxind] = max(avg);
                    p.thresh = avg(minind)+2*std(goodsnips(minind,:));
                    threshind = find(avg<p.thresh,1);
                    
                    p.win1del = minind-threshind;
                    p.win2del = maxind-threshind;
                    p.win1max = p.thresh;
                    p.win1min = avg(minind)-2*nanstd(goodsnips(minind,:));
                    p.win2max = avg(maxind)+2*nanstd(goodsnips(maxind,:));
                    p.win2min = max(0,avg(maxind)-2*nanstd(goodsnips(maxind,:)));
                    
                    if(any(structfun(@isempty, p)))
                        sortedspikes(chnind) = cell(1,length(chnind));
                        spikeshape(chnind) = cell(1,length(chnind));
                        continue;
                    end
                    
                    % Sort
                    detected = TwoWindowDiscrim(snips, p);
                    
                    spiketimes{1} = [spiketimes{1}, cross(detected)/fs];
                    spikeshape{chnind} = nanmean(goodsnips,2);
                    
                    if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                    
                    %% PCA if multiple spikes
                else
                    
                    if(verbose), tic; fprintf('  PCA...'); end
                    
                    % Sort
                    tempsnips = snips;
                    tempsnips(isnan(snips)) = 0;
                    
                    [coeff,score,~,~,explained,~] = pca(tempsnips');
                    variance = cumsum(explained/sum(explained));
                    r = find(variance > 0.90,1);
                    
                    Vr = score(:,1:r);
                    
                    idx = kmeans(Vr,nspikes,'Replicates',50,'MaxIter',1000);
                    
                    centroid = [];
                    
                    for i = 1:max(idx)
                        if(sum(idx==i) == 1)
                            centroid(i,:) = Vr(idx==i,:);
                        else
                            centroid(i,:) = mean(Vr(idx==i,:));
                        end
                        spiketimes{i} = [spiketimes{i}, cross(idx==i)/fs];
                        spikeshape{chnind(i)} = nanmean(snips(:,idx==i),2);
                    end
                    
                    if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                    
                end
                
                %% Loop through rest of data to find spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                start = 300;
                width = 300;
                while start < session_time
                    
                    if(verbose), tic; fprintf('  Loading data %d s...', start); end
                    
                    fin = start+width;
                    if(fin > session_time), fin = session_time; end
                    [wave, ~, ~, ~, ~, ~, ~, ~, ~, ~,~, session_time] = ...
                        LoadRippleData('path',path,'filename',filename,'flag',1,...
                        'wavechannel',allchns(c),'times',[start,fin]);
                    
                    % Filter data and get snippets
                    temp = HPF(wave,fs,916);
                    tempstim = SD.stim-start; tempstim(tempstim<0) = []; tempstim(tempstim>(length(wave)/fs)) = [];
                    [snips, cross] = GetSnips(temp, threshold, spikerange, 1000, tempstim,  0.0012, fs);
                    
                    if(isempty(snips))
                        start = fin;
                        if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                        continue;
                    end
                    
                    if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                    
                    % Sort
                    if(nspikes==1)
                        
                        if(verbose), tic; fprintf('  2 window discrim...'); end
                        
                        % 2 window
                        detected = TwoWindowDiscrim(snips, p);
                        
                        spiketimes{1} = [spiketimes{1}, cross(detected)/fs + start];
                        
                        if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                        
                    else
                        
                        if(verbose), tic; fprintf('  PCA...'); end
                        
                        % PCA
                        tempsnips = snips;
                        tempsnips(isnan(snips)) = 0;
                        
                        Vr = tempsnips'*coeff; Vr = Vr(:,1:r); Vr = Vr-mean(Vr);
                        
                        dist = pdist2(centroid,Vr); % distance to centroid
                        
                        [~,idx] = min(dist); % choose closest centroid
                        
                        for i = 1:max(idx)
                            spiketimes{i} = [spiketimes{i}, cross(idx==i)/fs + start];
                        end
                        
                        if(verbose), t = toc; fprintf(' %f seconds \n', t); end
                        
                    end
                    start = fin;
                end
                
                sortedspikes(chnind) = spiketimes;
                
                if(verbose), fprintf('Done\n'); end
                
            end
            
            % Fix channels if shifted
            if(str2num(days{d})<20190917 && max(chns)<64)
                chns = chns+32;
            end
            
            % Save
            if(verbose), tic; fprintf('  Saving...'); end
            save(fullfile(path,[filename,'_OfflineSort.mat']),'chns','codes','sortedspikes','spikeshape','-v7.3');
            if(verbose), t = toc; fprintf(' %f seconds \n\n', t); end
            
            
        catch
        end
        
    end
end

