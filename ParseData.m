%% Goes through and makes the session data struct for each experiment

%% Parse all data
clear; close all;

%% Load
% Path all folders are in
basepath = 'R:\Yun\Jafar\Ripple';

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
for d = 23%1:length(days)
    
    path = fullfile(basepath,days{d});
    files = dir(fullfile(path,'*.nev')); %find all ripple files
    files = extractfield(files,'name');
    
    for f = 1:length(files)
        
        %% Check to see if already parsed
        sdfile = fullfile(path,[files{f}(1:end-4),'_SD.mat']);
        if(exist(sdfile))
            continue;
        end
        
        disp([num2str(d),'. ',days{d},' - ',num2str(f),'. ',files{f}]);
        filename = files{f}(1:end-4);
        
        %% Initialize session data
        SD = struct;
        SD.Date = days{d};
        SD.Session = filename;
        
        try
            %% Load stim
            fprintf('Loading Stim Data...');
            tic;
            stimfile = fullfile(path,[files{f}(1:end-4),'_Stim.mat']);
            if(exist(stimfile))
                load(stimfile)
            end
            if(~exist(stimfile) || ~exist('stimchannel') || length(stimchannel)>1)
                [~, ~, ~, stim, stimchannel, ~, ~, ~, ~, ~,~]...
                    = LoadRippleData('path',path,'filename',filename,'flag',4);
                stim = stim.timestamp;
                if(length(stimchannel) > 1)
                    [~, ind] = max(sum(~isnan(stim),2));
                    stim = stim(ind,:);
                    stimchannel = stimchannel(ind);
                end
                save(stimfile,'stim','stimchannel');
            end
            SD.stim = stim;
            SD.stimchn = stimchannel;
            temp = SD.stim; temp = diff(temp); temp(temp<0.01) = [];
            SD.stimfreq = 1/mean(temp);
            t(1) = toc;
            fprintf('%1.2f seconds\n',t(1));
            
            %% Load spikes
            fprintf('Loading Spike Data...');
            tic;
            spikefile = fullfile(path,[filename,'_Spikes.mat']);
            if(exist(spikefile))
                load(spikefile)
            else
                [~, spikes, ~, ~, ~, ~, ~, ~, ~, ~,~]...
                    = LoadRippleData('path',path,'filename',filename,'flag',2);
                save(spikefile,'spikes');
            end
            t(2) = toc;
            fprintf('%1.2f seconds\n',t(2));
            
            %% Append analyzed session data
            fprintf('Parsing Data...');
            tic;
            % Evoked spikes
            [SD.chns, SD.codes, SD.ES, SD.ES_delay, SD.spiketimes] = extractES(spikes,stim);
            if(str2num(SD.Date)<20190917)
                SD.chns = SD.chns+32;
            end
            
            % Inhibition
            SD.IH = extractIH(SD.ES_delay,SD.spiketimes,stim);
            
            % Evoked spike trends
            [SD.ESslope, SD.ESpval] = extractTrend(SD.ES,stim,1);
            
            % Inhibition trends
            [SD.IHslope, SD.IHpval] = extractTrend(SD.IH,stim,0);
            
            % Channel distance
            SD.dist = [];
            [scol,srow,~] = GetChannelPosition(SD.stimchn);
            for c = 1:length(SD.chns)
                [col,row,~] = GetChannelPosition(SD.chns(c));
                SD.dist(c) = pdist([scol,srow;col,row])*0.4;
            end
            t(3) = toc;
            fprintf('%1.2f seconds\n',t(3));
            
            %% Save
            fprintf('Saving...');
            tic;
            save(sdfile,'SD');
            t(5) = toc;
            fprintf('%1.2f seconds\n',t(5));
            
            fprintf('Done in %1.2f seconds\n',sum(t));
            
        catch
            warning('Something happened');
        end
        
        % Clean some variables to prevent errors
        clearvars -except basepath days d path files f
    end
end





