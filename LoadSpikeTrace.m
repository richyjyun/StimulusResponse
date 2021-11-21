function [spikes] = LoadSpikeTrace(path,filename,chn)
% Summary
%   Returns spike traces of the denoted channel
%
% Inputs
%   path            path of ripple file
%   filename        ripple file name
%   chn             Utah channel
%
% Outputs
%   spikes          struct with timestamps, traces, sortcode

spikes = [];

%% Check file exists
fname = fullfile(path,[filename,'.nev']);
if(~exist(fname)), warning('No such file in directory'); return; end

%% Create a neuroshare file handle
[ns_RESULT, hFile] = ns_OpenFile(fname);
if(~strcmp(ns_RESULT,'ns_OK')), warning('Could not read file'); return; end

%% Find correct entity ID and load info
EntityIndices = find([hFile.Entity(:).ElectrodeID] == getRippleChn(chn));
EntityTypes = find(strcmp(extractfield(hFile.Entity,'EntityType'),'Segment'));
entityID = intersect(EntityIndices,EntityTypes);

if(isempty(entityID)), warning('No such entity'); return; end

[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID);

if(~strcmp(ns_RESULT,'ns_OK')), warning('Entity not correct'); return; end

%% Load all traces
maxload = 5000;
spikeEventTime_s = zeros(min(entityInfo.ItemCount,maxload),1);
unit_id = zeros(min(entityInfo.ItemCount,maxload),1);
for i = 1:min(entityInfo.ItemCount,maxload)
    [ns_RESULT, spikeEventTime_s(i), spikeWindowData(:,i), ~, unit_id(i)] = ns_GetSegmentData(hFile, entityID, i);
    if(~strcmp(ns_RESULT,'ns_OK')), warning('Could not load spike'); return; end
end

%% Assign spikes
spikes.timestamp = spikeEventTime_s;
spikes.snips = spikeWindowData;
spikes.sortcode = unit_id;

ns_RESULT = ns_CloseFile(hFile);
clear hFile

end

