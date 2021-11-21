%% Compiles all data from experiments into a single session list

%% Parse all data
clear; close all;

%% Load
% Path all folders are in
basepath = 'R:\Ripple';

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
SL = struct([]);
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
        
        % Remove bad days. 20190826 was after a ton of stim and lots of
        % noise/struggling. 
        if str2num(days{d})==20190826 && strcmp(files{f}(1:end-4),'20Hz_10uA_92')
            continue;
        end
        
        %% Check to see if SD exists
        sdfile = fullfile(path,[files{f}(1:end-4),'_SDNew.mat']);
        if(~exist(sdfile))
            continue;
        else
            load(sdfile);
        end
        
        if(max(SD.stim)-min(SD.stim) < 240)
            continue;
        end
        
        disp([num2str(d),'. ',days{d},' - ',num2str(f),'. ',files{f}]);
        filename = files{f}(1:end-4);
        
        SL = [SL,SD];
        
    end
end

slname = 'R:\Ripple\SL.mat';
save(slname,'SL','-v7.3');






