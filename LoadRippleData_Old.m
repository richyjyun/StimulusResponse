function [waveform, spikes, analog, stim, stimchannel, wavechannel, analogchannel, spikechannel, fs, times,lfp, session_time]...
    = LoadRippleData_Old(varargin)
% Summary
%   Loads data from Neuroshare file format for Ripple data. Need Neuroshare
%   functions to operate.
%
% Inputs
%   path            path of ripple file
%   filename        ripple file name
%   flag            what data to return
%                   0 all data (default)
%                   1 raw waveform data 
%                   2 spike data
%                   3 analog input data
%                   4 stim data
%                   5 lfp data
%   wavechannel     returning specific channels (Utah channels) if denoted, otherwise return all
%   analogchannel   returning specific channels if denoted, otherwise return all
%   spikechannel    returning specific channels if denoted, otherwise return all
%   lfpchannel      returning specific channels if denoted, otherwise return all
%   times           time epoch to return in seconds. empty to return all
%
% Outputs
%   waveform        waveform data matrix
%   spikes          struct with timestamps, channel, and sortcode of spikes
%   analog          analog input data matrix
%   fs              sampling rate
%   channel         channel returned, NaN if all channels
%   times           time epoch returned

%% Set up name/value pairs

p = inputParser;

addParameter(p,'path',cd);
addParameter(p,'filename',[]);
addParameter(p,'flag',0);
addParameter(p,'wavechannel',nan);
addParameter(p,'analogchannel',nan);
addParameter(p,'spikechannel',nan);
addParameter(p,'lfpchannel',nan);
addParameter(p,'times',[]);

parse(p,varargin{:});

path = p.Results.path;
filename = p.Results.filename;
flag = p.Results.flag;
wavechannel = p.Results.wavechannel;
analogchannel = p.Results.analogchannel;
spikechannel = p.Results.spikechannel;
lfpchannel = p.Results.lfpchannel;
times = p.Results.times;

waveform = [];
spikes = [];
analog = [];
lfp = [];
stim = [];
stimchannel = [];

%% Input error checking
if(isempty(filename)), error('No filename given'); end

if(size(path,1)>1)
    path = path';
end
fname = fullfile(path,[filename,'.nev']);

if(~exist(fname)), error('No such file in directory'); end
if(any(flag<0) || any(flag>5)), warning('Flag outside range, returning all data'); flag = 0; end

%% Create a neuroshare file handle
[ns_RESULT, hFile] = ns_OpenFile(fname);
if(~strcmp(ns_RESULT,'ns_OK')), error('Could not read file'); end

%% Read high-level information about the recording
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
session_time = nsFileInfo.TimeSpan;
if(~strcmp(ns_RESULT,'ns_OK')), error('Could not read file info'); end

if(isempty(times))
    times=[0,nsFileInfo.TimeSpan];
else
    if(times(1) > times(2) || length(times)~=2), error('Time denoted not proper'); end
    if(times(1)>nsFileInfo.TimeSpan), error('Time epoch is outside of time span'); end
    if(times(2)>nsFileInfo.TimeSpan), times(2)=nsFileInfo.TimeSpan; warning('End of epoch changed to end of time span'); end
end

fs = 1/nsFileInfo.TimeStampResolution;

% Preallocate memory for nsEntityInfo
nsEntityInfo(nsFileInfo.EntityCount,1).EntityLabel = '';
nsEntityInfo(nsFileInfo.EntityCount,1).EntityType  = 0;
nsEntityInfo(nsFileInfo.EntityCount,1).ItemCount   = 0;

% The entity inforamtion is read using ns_GetEntityInfo()
for i = 1:nsFileInfo.EntityCount
    [~, nsEntityInfo(i,1)] = ns_GetEntityInfo(hFile, i);
end

% If flag denotes all data or waveform or analog
if(any(ismember(flag,[0,1,3,4,5])))
    
    %% Determine analog entities and labels
    %Get stream entities
    EntityIDs  = find([nsEntityInfo.EntityType]==2);
    SegmentEntityIDs = find([nsEntityInfo.EntityType]==3);

    % Find labels
    EntityLabels  = {nsEntityInfo(EntityIDs).EntityLabel};
    SegmentEntityLabels = {nsEntityInfo(SegmentEntityIDs).EntityLabel};

    % Organize item counts by entity type
    ItemCounts  = [nsEntityInfo(EntityIDs).ItemCount];
    SegmentItemCounts = [nsEntityInfo(SegmentEntityIDs).ItemCount];

    %% Split between analog input and other
    analogind = contains(EntityLabels,'analog');
    waveind = contains(EntityLabels,'raw');
    lfpind = contains(EntityLabels,'lfp');
    
    AnalogEntityIDs = EntityIDs(analogind);
    AnalogEntityLabels = EntityLabels(analogind);
    AnalogItemCounts = ItemCounts(analogind);
    
    WaveEntityIDs = EntityIDs(waveind);
    WaveEntityLabels = EntityLabels(waveind);
    WaveItemCounts = ItemCounts(waveind);
    
    LFPEntityIDs = EntityIDs(lfpind);
    LFPEntityLabels = EntityLabels(lfpind);
    LFPItemCounts = ItemCounts(lfpind);
    
    ReadStart = round(times(1)*fs); if(ReadStart <= 0), times(1) = 0; ReadStart = 1; end
    SampleCount = round((times(2)-times(1))*fs);
    
    %% Get waveform data
    if(any(ismember(flag,[0,1])))
        
        WaveEntityLabels = cellfun(@(x) str2num(x(5:end)),WaveEntityLabels);
        WaveEntityLabels = getUtahChn(WaveEntityLabels);
        
        WaveEntityCount = length(WaveEntityIDs);
        WaveSampleReadCounts = min(WaveItemCounts, SampleCount);
        fprintf('Retrieving wave data for channel: ');
        
        % Get all channels if not denoted, otherwise find correct index
        if(isnan(wavechannel)) 
            wavechannel = 1:WaveEntityCount; 
        else
            wavechannel = find(ismember(WaveEntityLabels,wavechannel));
        end
        
        waveform = zeros(max(WaveSampleReadCounts), length(wavechannel));

        % Read data
        for i = wavechannel
            display(num2str(WaveEntityLabels(i)));
            [~, ~, waveform(1:WaveSampleReadCounts(i),find(wavechannel==i))] = ...
                ns_GetAnalogData(hFile, WaveEntityIDs(i), ...
                ReadStart, WaveSampleReadCounts(i),'');
        end        
    end
    
    %% Get analog input data
    if(any(ismember(flag,[0,3])))
        AnalogEntityCount = length(AnalogEntityIDs);
        AnalogSampleReadCounts = min(AnalogItemCounts, SampleCount);
        fprintf('Retrieving analog data for channel:\n');
        
        % Get all channels if not denoted
        if(isnan(analogchannel)), analogchannel = 1:AnalogEntityCount; end
        
        analog = zeros(max(AnalogSampleReadCounts), length(analogchannel));
        
        for i = analogchannel
            display(num2str(i));
            % Read data 
            [~, ~, analog(1:AnalogSampleReadCounts(i),find(analogchannel==i))] = ...
                ns_GetAnalogData(hFile, AnalogEntityIDs(i), ...
                ReadStart, AnalogSampleReadCounts(i),'');
        end        
    end
    
    %% Get LFP data
    if(any(ismember(flag,[0,5])))
        
        % LFPs are read in at 1kHz
        ReadStart = round(times(1)*1000); if(ReadStart <= 0), times(1) = 0; ReadStart = 1; end
        SampleCount = round((times(2)-times(1))*1000);
        
        LFPEntityLabels = cellfun(@(x) str2num(x(5:end)),LFPEntityLabels);
        LFPEntityLabels = getUtahChn(LFPEntityLabels);
        
        LFPEntityCount = length(LFPEntityIDs);
        LFPSampleReadCounts = min(LFPItemCounts, SampleCount);
        fprintf('Retrieving wave data for channel: ');
        
        % Get all channels if not denoted, otherwise find correct index
        if(isnan(lfpchannel)) 
            lfpchannel = 1:LFPEntityCount; 
        else
            lfpchannel = find(ismember(LFPEntityLabels,lfpchannel));
        end
        
        lfp = zeros(max(LFPSampleReadCounts), length(lfpchannel));

        % Read data
        for i = lfpchannel
            display(num2str(WaveEntityLabels(i)));
            [~, ~, lfp(1:LFPSampleReadCounts(i),find(lfpchannel==i))] = ...
                ns_GetAnalogData(hFile, LFPEntityIDs(i), ...
                ReadStart, LFPSampleReadCounts(i),'');
        end        
    end
    
    
    %% Get stim data
    if(any(ismember(flag,[0,4])))
        IDs = extractfield(hFile.Entity,'ElectrodeID');
        Types = extractfield(hFile.Entity,'EntityType');
        IDs = IDs(strcmp(Types,'Segment'));
        
        ind = IDs>5120; % Electrode assignment for stim
        
        SegmentEntityCount   = length(SegmentEntityIDs(ind));
        SegmentTimeStamps    = nan(SegmentEntityCount, max(SegmentItemCounts(ind)));
        fprintf('Retrieving stim\n')
        
        inds = find(ind);
        
        stimchannel = getUtahChn(IDs(inds)-5120);
        
        for i = 1:length(inds)
            display(num2str(getUtahChn(IDs(inds(i))-5120)));

            % Get stim
            for j = 1:SegmentItemCounts(inds(i))
                [~,SegmentTimeStamps(i,j), ~, ~, ~] = ...
                    ns_GetSegmentData(hFile, SegmentEntityIDs(inds(i)), j);
            end
        end

        remove = all(isnan(SegmentTimeStamps),2);
        
        stim.labels = SegmentEntityLabels(inds);
        stim.timestamp = SegmentTimeStamps;
        
    end
    
end

if(any(ismember(flag,[0,2])))
        fprintf('Retrieving spiking data\n');
    try
        %% Get timestamps for spikes
        TS = hFile.FileInfo(1).MemoryMap.Data.TimeStamp;
        CHN = hFile.FileInfo(1).MemoryMap.Data.PacketID;
        CHN = getUtahChn(CHN);
        SC = hFile.FileInfo(1).MemoryMap.Data.Class;
        
        bad = TS < times(1)*fs | TS > times(2)*fs | CHN == spikechannel | SC == 0;
        
        spikes.timestamp = double(TS(~bad))./fs;
        spikes.chan = double(CHN(~bad));
        spikes.sortcode = double(SC(~bad));
    catch
    end
end

%% Close the file handle
ns_RESULT = ns_CloseFile(hFile);
if( ~strcmp(ns_RESULT, 'ns_OK')), disp(['ERROR: ns_CloseFile() returned ' ns_RESULT]); end
clear hFile

end
