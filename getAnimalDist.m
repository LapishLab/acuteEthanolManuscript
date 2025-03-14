% Load the data
load('allAnimalsDistData.mat', 'allBinnedDistance');

% Identify ET and EN animals/recordings
etCols = endsWith(allBinnedDistance(211, :), 'et');
enCols = endsWith(allBinnedDistance(211, :), 'en');

% Extract and convert to array
distEt = cell2mat(allBinnedDistance(1:210, etCols));
distEn = cell2mat(allBinnedDistance(1:210, enCols));

% Add animal numbers as the last row for ID
distEt = [distEt; cellfun(@(x) str2double(extractBefore(x, 'et')), allBinnedDistance(211, etCols))];
distEn = [distEn; cellfun(@(x) str2double(extractBefore(x, 'en')), allBinnedDistance(211, enCols))];

