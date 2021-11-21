%% Editing SL if certain algorithms were changed

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

%% Go through all the data
for d = 1:length(days)
    
    path = fullfile(basepath,days{d});
    files = dir(fullfile(path,'*.nev')); %find all ripple files
    files = extractfield(files,'name');
    
    for f = 1:length(files)
        
        %% Check to see if relevant file
        if(~(contains(files{f},'5Hz') || contains(files{f},'10Hz') || contains(files{f},'20Hz') ||...
                contains(files{f},'Gamm')))
            continue;
        end
        
        %% Check to see if SD exists
        sdfile = fullfile(path,[files{f}(1:end-4),'_SD.mat']);
        if(~exist(sdfile))
            continue;
        else
            load(sdfile);
        end
        
        disp([num2str(d),'. ',days{d},' - ',num2str(f),'. ',files{f}]);
        filename = files{f}(1:end-4);
        
        t = [];
%         try
            %% Load spikes 
            fprintf('Loading Spike Data...');
            tic;
            
            % Online spike sorting
            spikefile = fullfile(path,[filename,'_Spikes.mat']);
            if(exist(spikefile))
                load(spikefile)
            else
                [~, spikes, ~, ~, ~, ~, ~, ~, ~, ~,~]...
                    = LoadRippleData('path',path,'filename',filename,'flag',2);
                save(spikefile,'spikes');
            end
            
%             % Offline spike sorting. Reformat for extractES
%             spikefile = fullfile(path,[filename,'_OfflineSort.mat']);
%             load(spikefile);
%             spikes = struct; spikes.chan = []; spikes.sortcode = []; spikes.timestamp = [];
%             for i = 1:length(sortedspikes)
%                 spikes.timestamp = [spikes.timestamp,sortedspikes{i}];
%                 spikes.chan = [spikes.chan, chns(i)*ones(1,length(sortedspikes{i}))];
%                 spikes.sortcode = [spikes.sortcode, codes(i)*ones(1,length(sortedspikes{i}))];
%             end
%             spikes.timestamp = spikes.timestamp';
%             spikes.chan = spikes.chan';
%             spikes.sortcode = spikes.sortcode';
            
            t(2) = toc;
            fprintf('%1.2f seconds\n',t(2));
            
            %% Append analyzed session data
            fprintf('Parsing Data...');
            tic;
            
            % Fix stim
            temp = SD.stim; temp = diff(temp); % remove stim that are too close
            SD.stim(temp<0.01) = [];
            
            % Fix evoked spikes
            [SD.chns, SD.codes, SD.ES, SD.ES_delay, SD.spiketimes] = extractES_New(spikes,SD.stim);
            if(str2num(SD.Date)<20190917 && max(SD.chns)<64)
                SD.chns = SD.chns+32;
            end
            
            % Inhibition
%             SD.IH = extractIH(SD.ES_delay,SD.spiketimes,SD.stim);
            
            % Evoked spike trends
            [SD.ESslope, SD.ESpval] = extractTrend(SD.ES,SD.stim,1);
            
            % Inhibition trends
%             [SD.IHslope, SD.IHpval] = extractTrend(SD.IH,SD.stim,0);
            
            % Channel distance
            SD.dist = [];
            [scol,srow,~] = GetWadeChannelPosition(SD.stimchn);
            for c = 1:length(SD.chns)
                [col,row,~] = GetWadeChannelPosition(SD.chns(c));
                SD.dist(c) = pdist([scol,srow;col,row])*0.4;
            end
            
            t(3) = toc;
            fprintf('%1.2f seconds\n',t(3));
            
            %% Save
            fprintf('Saving...');
            tic;
            savefile = fullfile(path,[files{f}(1:end-4),'_SDNew.mat']);
            save(savefile,'SD');
            t(end+1) = toc;
            fprintf('%1.2f seconds\n',t(end));
            
            fprintf('Done in %1.2f seconds\n',sum(t));
            
%         catch
%             warning('Something happened');
%         end
        
        % Clean some variables to prevent errors
        clearvars -except basepath days d path files f
    end
end

MakeSL();


