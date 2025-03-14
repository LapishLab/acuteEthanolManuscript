%% ratAcuteEthanol.m
% This code generates the data/graphs for Figure 4 pertaining to the
% zscored mean firing, percent correlated, and zscored firing for
% positively and negatively correlated populations for Ensure only and
% Ensure + ETOH

clear all
close all
% only pull .mat files that have en or et in them
% en = ensure only recording
% et = ensure + ethanol recording
fileList = dir('*.mat');
pattern = '.*(en|et).*\.mat';
matchingFiles = arrayfun(@(x) ~isempty(regexp(x.name, pattern, 'once')), fileList);
fileList = fileList(matchingFiles);

%% Pull out timestamps, spikes, and bin
% This pulls out spikes separately as Good and MUA, then at the end
% of the script they get concatenated for a giant block of 
% good + MUA units for ensure only and ensure + ethanol
% Modified this from binnedArrays.m to only look at data after the
% bottle is added
% Figure out max neuron numbers
maxNeuronsGood = 0;
maxNeuronsMUA = 0;

for iFile = 1:length(fileList)
    file = load(fileList(iFile).name);
    numNeuronsGood = size(file.spkGood, 1);
    numNeuronsMUA = size(file.spkMUA, 1);
    maxNeuronsGood = max(maxNeuronsGood, numNeuronsGood);
    maxNeuronsMUA = max(maxNeuronsMUA, numNeuronsMUA);
end

% Bin size variables here
binUnder1800s = 30;
binOver1800s = 180;

% Initialize arrays for different conditions and bins
binnedSpikesEnGood = NaN(maxNeuronsGood * length(fileList), binUnder1800s + binOver1800s);
binnedSpikesEtGood = NaN(maxNeuronsGood * length(fileList), binUnder1800s + binOver1800s);
binnedSpikesEnMUA = NaN(maxNeuronsMUA * length(fileList), binUnder1800s + binOver1800s);
binnedSpikesEtMUA = NaN(maxNeuronsMUA * length(fileList), binUnder1800s + binOver1800s);

animalNumberEtGood = zeros(maxNeuronsGood * length(fileList), 1); 
animalNumberEtMUA = zeros(maxNeuronsMUA * length(fileList), 1); 
animalNumberEnGood = zeros(maxNeuronsGood * length(fileList), 1); 
animalNumberEnMUA = zeros(maxNeuronsMUA * length(fileList), 1); 

for iFile = 1:length(fileList)
    file = load(fileList(iFile).name);
    
    animalNumber = str2double(regexp(fileList(iFile).name, '(\d+)(et|en)\.mat', 'tokens', 'once'));

    numNeuronsGood = size(file.spkGood, 1);
    numNeuronsMUA = size(file.spkMUA, 1);
% Extrapolate recording time from earliest/latest spike in each file
% There are no explicit recording times in the .mat files (that i see)
% So, using earliest and latest time stamps to extrapolate. They are all 
% pretty close if you look at the recordingTimes (~12,000s)
% This gets used to determine the recording time for each file and
% the earliest time stamps also get used to determine time zero
% from which the 1800s (i.e., 30 min) baseline period is determined

    spkGoodTimestamps = cellfun(@(x) x, file.spkGood, 'UniformOutput', false);
    spkMUATimestamps = cellfun(@(x) x, file.spkMUA, 'UniformOutput', false);
    allTimestamps = [];
    for i = 1:length(spkGoodTimestamps)
        allTimestamps = [allTimestamps; spkGoodTimestamps{i}];
    end
    for i = 1:length(spkMUATimestamps)
        allTimestamps = [allTimestamps; spkMUATimestamps{i}];
    end
    
    % Find the earliest and latest
    earliestTimestamp = min(allTimestamps);
    latestTimestamp = max(allTimestamps);
    
    % Extrapolate the recording time for the entire file
    recordingTime = latestTimestamp - earliestTimestamp;

    % Then, check for "en" or "et" in the file name
    % and save them in separate arrays as appropriate
    if contains(fileList(iFile).name, 'en')
        % Iterate through spkGood neurons
        for iNeuron = 1:numNeuronsGood
            spikeTimes = file.spkGood{iNeuron}; % get spike times
            index = (iFile - 1) * maxNeuronsGood + iNeuron; % figure out index
            
            % Divide spike times into two segments based on 1800s (30 min) threshold
            spikeUnder1800s = spikeTimes(spikeTimes < 1800 + earliestTimestamp); % using the earliestTimestamp here as time zero for the whole file
            spikeOver1800s = spikeTimes(spikeTimes >= 1800 + earliestTimestamp) - (1800 + earliestTimestamp); % get everyhting after 1800s
            [counts, ~] = histcounts(spikeUnder1800s, binUnder1800s); %bin the spikes pre 1800s
            binWidth = 60; % 60 because 1800 (time period) /30 (bin #) is 60
            firingRate = counts / binWidth; % determine FRs for each bin from the time period
            binnedSpikesEnGood(index, 1:binUnder1800s) = firingRate; % insert into the array
            [counts, ~] = histcounts(spikeOver1800s, binOver1800s); % same but for the rest of the spikes after 1800s
            binWidth = (recordingTime - 1800) / binOver1800s; % this will vary per file because of the slightly different recording end times
            firingRate = counts / binWidth;
            binnedSpikesEnGood(index, binUnder1800s+1:end) = firingRate;% insert into array after the 30 bins
            animalNumberEnGood(index) = animalNumber(1);
        end
        
        % Same thing but for spkMUA now
        for iNeuron = 1:numNeuronsMUA
            spikeTimes = file.spkMUA{iNeuron};
            index = (iFile - 1) * maxNeuronsMUA + iNeuron;
            spikeUnder1800s = spikeTimes(spikeTimes < 1800 + earliestTimestamp);
            spikeOver1800s = spikeTimes(spikeTimes >= 1800 + earliestTimestamp) - (1800 + earliestTimestamp);
            [counts, ~] = histcounts(spikeUnder1800s, binUnder1800s);
            binWidth = 60;
            firingRate = counts / binWidth;
            binnedSpikesEnMUA(index, 1:binUnder1800s) = firingRate;
            [counts, ~] = histcounts(spikeOver1800s, binOver1800s);
            binWidth = (recordingTime - 1800) / binOver1800s;
            firingRate = counts / binWidth;
            binnedSpikesEnMUA(index, binUnder1800s+1:end) = firingRate;
            animalNumberEnMUA(index) = animalNumber(1);
        end
    elseif contains(fileList(iFile).name, 'et')
        % Same as above but for the ethanol files now
        for iNeuron = 1:numNeuronsGood
            spikeTimes = file.spkGood{iNeuron};
            index = (iFile - 1) * maxNeuronsGood + iNeuron;
            spikeUnder1800s = spikeTimes(spikeTimes < 1800 + earliestTimestamp);
            spikeOver1800s = spikeTimes(spikeTimes >= 1800 + earliestTimestamp) - (1800 + earliestTimestamp); 
            [counts, ~] = histcounts(spikeUnder1800s, binUnder1800s);
            binWidth = 60;
            firingRate = counts / binWidth;
            binnedSpikesEtGood(index, 1:binUnder1800s) = firingRate;
            [counts, ~] = histcounts(spikeOver1800s, binOver1800s);
            binWidth = (recordingTime - 1800) / binOver1800s;
            firingRate = counts / binWidth;
            binnedSpikesEtGood(index, binUnder1800s+1:end) = firingRate;
            animalNumberEtGood(index) = animalNumber(1);
        end
        
        for iNeuron = 1:numNeuronsMUA
            spikeTimes = file.spkMUA{iNeuron};
            index = (iFile - 1) * maxNeuronsMUA + iNeuron;
            spikeUnder1800s = spikeTimes(spikeTimes < 1800 + earliestTimestamp);
            spikeOver1800s = spikeTimes(spikeTimes >= 1800 + earliestTimestamp) - (1800 + earliestTimestamp); 
            [counts, ~] = histcounts(spikeUnder1800s, binUnder1800s);
            binWidth = 60;
            firingRate = counts / binWidth;
            binnedSpikesEtMUA(index, 1:binUnder1800s) = firingRate;
            [counts, ~] = histcounts(spikeOver1800s, binOver1800s);
            binWidth = (recordingTime - 1800) / binOver1800s;
            firingRate = counts / binWidth;
            binnedSpikesEtMUA(index, binUnder1800s+1:end) = firingRate;
            animalNumberEtMUA(index) = animalNumber(1);
        end
    end
end
binnedSpikesEnGood(:, end+1) = animalNumberEnGood;
binnedSpikesEnMUA(:, end+1) = animalNumberEnMUA;
binnedSpikesEtGood(:, end+1) = animalNumberEtGood;
binnedSpikesEtMUA(:, end+1) = animalNumberEtMUA;

%% Remove the rows with NaN values
binnedSpikesEtGood = binnedSpikesEtGood(~any(isnan(binnedSpikesEtGood), 2), :);
binnedSpikesEtMUA = binnedSpikesEtMUA(~any(isnan(binnedSpikesEtMUA), 2), :);
binnedSpikesEnGood = binnedSpikesEnGood(~any(isnan(binnedSpikesEnGood), 2), :);
binnedSpikesEnMUA = binnedSpikesEnMUA(~any(isnan(binnedSpikesEnMUA), 2), :);
concatSpikesEt = [binnedSpikesEtGood; binnedSpikesEtMUA];
concatSpikesEn = [binnedSpikesEnGood; binnedSpikesEnMUA];
% delete row/neuron 98 which seems to be noise based on imagesc/zscore
concatSpikesEt(98,:) = [];
% Sort based on animal numbers
animalNumSpikesEt = sortrows(concatSpikesEt, size(concatSpikesEt, 2));
animalNumSpikesEn = sortrows(concatSpikesEn, size(concatSpikesEn, 2));
%% Zscore
zSpikesEt = zscore(animalNumSpikesEt(:,1:210)');
baseZEt = zSpikesEt(1:30,:);
baseZEtMean = mean(baseZEt);
zSpikesEt = zSpikesEt - baseZEtMean;

zSpikesEn = zscore(animalNumSpikesEn(:,1:210)');
baseZEn = zSpikesEn(1:30,:);
baseZEnMean = mean(baseZEn);
zSpikesEn = zSpikesEn - baseZEnMean;

animalNumSpikesEt(:,1:210) = zSpikesEt';
animalNumSpikesEn(:,1:210) = zSpikesEn';
%% Plot mean zscored FRs with filled in standard error bars
meanEt = mean(concatSpikesEt);
meanEn = mean(concatSpikesEn);

steEt = std(concatSpikesEt) / sqrt(size(concatSpikesEt, 1));
steEn = std(concatSpikesEn) / sqrt(size(concatSpikesEn, 1));

stdEt = std(concatSpikesEt);
stdEn = std(concatSpikesEn);

figure;
plot(meanEn, 'r', 'LineWidth', 0.5);
hold on;
plot(meanEt, 'b', 'LineWidth', 0.5);
%fill_between function is at the bottom - used multiple times in this
%script
fill_between(1:length(meanEt), meanEt - steEt, meanEt + steEt, 'b', 0.3);
fill_between(1:length(meanEn), meanEn - steEn, meanEn + steEn, 'r', 0.3);
xlabel('Time');
ylabel('Zscored mean FR');
xlim([0, 210]);
legend('Ensure Only','Ensure + Ethanol');
xline(30,'g','Bottle Added')

%% Load intake data and prep/sort
interpolateIntakes; %calls the interpolateIntakes script. Requires the ephysensurenoncumul.xlsx and ephysetohnoncumul.xlsx sheet for this
newColumn = [16, 18, 19, 37, 38, 17]; %their corresponding subject numbers for Etoh
numColumns = size(itpEtIntake, 2) + 1;
itpEtIntake(:, numColumns) = newColumn;

newColumn = [16, 17, 18, 19, 37, 38]; %corresponding subject numbers for ensure
numColumns = size(itpEnIntake, 2) + 1;
itpEnIntake(:, numColumns) = newColumn;
% Sort intakes based on animal numbers
animalNumEtIntake = sortrows(itpEtIntake, size(itpEtIntake,2));
animalNumEnIntake = sortrows(itpEnIntake, size(itpEnIntake,2));

%% Correlate each neuron with each animals intake
% For ethanol
numNeurons = size(animalNumSpikesEt, 1);
numAnimals = size(animalNumEtIntake, 1);

correlEt = zeros(numNeurons, 2); 

for iNeuron = 1:numNeurons
    neuronData = animalNumSpikesEt(iNeuron, 31:210);
    
    for iAnimal = 1:numAnimals
        animalID = animalNumSpikesEt(iNeuron, 211);
        
        % Find the corresponding intake data for the current animal
        intakeData = animalNumEtIntake(animalNumEtIntake(:, 181) == animalID, 1:180);
        
        % Check if the number of rows matches
        if size(neuronData, 2) == size(intakeData, 2)
            % Calculate correlation coefficient and p-value
            [corrCoeff, pValue] = corrcoef(neuronData, intakeData);
            
            % Store correlation coefficient and p-value in the array
            correlEt(iNeuron, 1) = corrCoeff(1, 2); % Correlation coefficient
            correlEt(iNeuron, 2) = pValue(1, 2);    % P-value
        end
    end
end

% Same but for ensure
numNeurons = size(animalNumSpikesEn, 1);
numAnimals = size(animalNumEnIntake, 1);

correlEn = zeros(numNeurons, 2); 
for iNeuron = 1:numNeurons
    neuronData = animalNumSpikesEn(iNeuron, 31:210);
    
    for iAnimal = 1:numAnimals
        animalID = animalNumSpikesEn(iNeuron, 211);
        
        % Find the corresponding intake data for the current animal
        intakeData = animalNumEnIntake(animalNumEnIntake(:, 181) == animalID, 1:180);
        
        % Check if the number of rows matches
        if size(neuronData, 2) == size(intakeData, 2)
            % Calculate correlation coefficient and p-value
            [corrCoeff, pValue] = corrcoef(neuronData, intakeData);
            
            % Store correlation coefficient and p-value in the array
            correlEn(iNeuron, 1) = corrCoeff(1, 2); % Correlation coefficient
            correlEn(iNeuron, 2) = pValue(1, 2);    % P-value
        end
    end
end
%% Find signficant correlations
%Find and index signficant positive and negative correlations for Ensure
alpha = 0.05;
sigPEn = correlEn(:, 2) < alpha;
sigR2En = correlEn(sigPEn, 1);

posIndEn = sigR2En > 0;
negIndEn = sigR2En < 0;
totalInd = length(correlEn);
%Determine percents from the total
percentPosInd = sum(posIndEn) / totalInd * 100;
percentNegInd = sum(negIndEn) / totalInd * 100;
percentNonResp = 100 - (percentPosInd + percentNegInd);
percentPosNegEn = [percentPosInd, percentNegInd, percentNonResp];
%Graph in a pie chart
figure;
subplot(1,2,1);
pie(percentPosNegEn);
legend({'Positive', 'Negative', 'No Correl'}, 'Location', 'Best');
title('Correlations for Ensure Only')

%Same as above but with Ethanol
sigPEt = correlEt(:, 2) < alpha;
sigR2Et = correlEt(sigPEt, 1);
posIndEt = sigR2Et > 0;
negIndEt = sigR2Et < 0;
totalInd = length(correlEt);

percentPosInd = sum(posIndEt) / totalInd * 100;
percentNegInd = sum(negIndEt) / totalInd * 100;
percentNonResp = 100 - (percentPosInd + percentNegInd);
percentPosNegEt = [percentPosInd, percentNegInd, percentNonResp];
subplot(1,2,2);
pie(percentPosNegEt);
legend({'Positive', 'Negative', 'No Correl'}, 'Location', 'Best');
title('Correlations for Ethanol')
%% Mean FR graphs depending on correlation
% Find indices with positive and negative r2 values for EtOH
positiveR2IndicesEt = find(correlEt(:, 1) > 0);
negativeR2IndicesEt = find(correlEt(:, 1) < 0);

significantIndicesEt = find(correlEt(:, 2) < 0.05);
positiveR2IndicesEt = intersect(positiveR2IndicesEt, significantIndicesEt);
negativeR2IndicesEt = intersect(negativeR2IndicesEt, significantIndicesEt);
positiveSpikesDataEt = animalNumSpikesEt(positiveR2IndicesEt, 1:210);
negativeSpikesDataEt = animalNumSpikesEt(negativeR2IndicesEt, 1:210);
%Calculate mean and STE
meanPositiveEt = mean(positiveSpikesDataEt);
stePositiveEt = std(positiveSpikesDataEt) / sqrt(size(positiveSpikesDataEt, 1));
meanNegativeEt = mean(negativeSpikesDataEt);
steNegativeEt = std(negativeSpikesDataEt) / sqrt(size(negativeSpikesDataEt, 1));


% Same as above but with Ensure
positiveR2IndicesEn = find(correlEn(:, 1) > 0);
negativeR2IndicesEn = find(correlEn(:, 1) < 0);

significantIndicesEn = find(correlEn(:, 2) < 0.05);
positiveR2IndicesEn = intersect(positiveR2IndicesEn, significantIndicesEn);
negativeR2IndicesEn = intersect(negativeR2IndicesEn, significantIndicesEn);
positiveSpikesDataEn = animalNumSpikesEn(positiveR2IndicesEn, 1:210);
negativeSpikesDataEn = animalNumSpikesEn(negativeR2IndicesEn, 1:210);

meanPositiveEn = mean(positiveSpikesDataEn);
stePositiveEn = std(positiveSpikesDataEn) / sqrt(size(positiveSpikesDataEn, 1));
meanNegativeEn = mean(negativeSpikesDataEn);
steNegativeEn = std(negativeSpikesDataEn) / sqrt(size(negativeSpikesDataEn, 1));
%% FDR post-hoc t-test
%MATLAB is bad at ANOVAs, so the rmANOVAs were performed with the exported
%data in Graphpad Prism. FDR post-hocs were then completed here.

%Calculate sig t-tests for Ethanol
pvalsEt = zeros(1, size(positiveSpikesDataEt, 2));

% Loop through each column
for i = 1:size(positiveSpikesDataEt, 2)
    % Extract the current column from 'pos' and 'neg'
    posSpikes = positiveSpikesDataEt(:, i);
    negSpikes = negativeSpikesDataEt(:, i);
    
    % Perform two-sample t-test
    [~, p_value] = ttest2(posSpikes, negSpikes);
    
    % Store the p-value for the current column
    pvalsEt(i) = p_value;
end
% Same as above but for Ensure
pvalsEn = zeros(1, size(positiveSpikesDataEn, 2));

% Loop through each column
for i = 1:size(positiveSpikesDataEn, 2)
    % Extract the current column from 'pos' and 'neg'
    posSpikes = positiveSpikesDataEn(:, i);
    negSpikes = negativeSpikesDataEn(:, i);
    
    % Perform two-sample t-test
    [~, p_value] = ttest2(posSpikes, negSpikes);
    
    % Store the p-value for the current column
    pvalsEn(i) = p_value;
end

%Calculate sig FDR post-hocs for EtOH and Ensure by calling the fdr_bh script
[hEt, crit_pEt, adj_ci_cvrgEt, adj_pEt]=fdr_bh(pvalsEt,0.05,'pdep','no');
[hEn, crit_pEn, adj_ci_cvrgEn, adj_pEn]=fdr_bh(pvalsEn,0.05,'pdep','no');

%% Plot the positively and negatively correlated neurons for EtOH and Ensure only
%with FDR sig. post-hocs signified by asteriks

% For ETOH
figure;
plot(meanPositiveEt, 'y', 'LineWidth', 2);
hold on;
plot(meanNegativeEt, 'r', 'LineWidth', 2);
fill_between(1:length(meanPositiveEt), meanPositiveEt - stePositiveEt, meanPositiveEt + stePositiveEt, 'y', 0.3);
fill_between(1:length(meanNegativeEt), meanNegativeEt - steNegativeEt, meanNegativeEt + steNegativeEt, 'r', 0.3);

ylabel('Zscored FR')
xlabel('Time')
title('Ensure + EtOH')
legend('Pos. corr.', 'Neg. corr.')
xlim([0, 210])
xFDR = find(hEt); 
yFDR = 1; 
hold on;
plot(xFDR, yFDR * ones(size(xFDR)), 'w*', 'MarkerSize', 10);

% for Ensure
figure
plot(meanPositiveEn, 'y', 'LineWidth', 2);
hold on;
plot(meanNegativeEn, 'r', 'LineWidth', 2);
fill_between(1:length(meanPositiveEn), meanPositiveEn - stePositiveEn, meanPositiveEn + stePositiveEn, 'y', 0.3);
fill_between(1:length(meanNegativeEn), meanNegativeEn - steNegativeEn, meanNegativeEn + steNegativeEn, 'r', 0.3);
ylabel('Zscored FR')
xlabel('Time')
title('Ensure only')
legend('Pos. corr.', 'Neg. corr.')
xlim([0, 210])
xFDR = find(hEn); 
yFDR = 2; 
hold on;
plot(xFDR, yFDR * ones(size(xFDR)), 'w*', 'MarkerSize', 10);
%% Graph mean distances for groups 
getAnimalDist; %run script which takes distances and attaches animal IDs to them in similar fashion as we have for neurons
distEt = distEt';%flip them so animals are in rows instead of columns
distEn = distEn';

% Calculate the mean for columns 1:210 across all animals (rows)
meanDistEt = mean(distEt(:, 1:210), 1);
meanDistEn = mean(distEn(:, 1:210), 1);

% Plot the means
figure;
plot(meanDistEt, 'r', 'LineWidth', 2); hold on;
plot(meanDistEn, 'y', 'LineWidth', 2);
hold off;

% Add labels and legend
xlabel('Time (minutes)');
ylabel('Mean Distance (cm)');
legend('Ensure + EtOH', 'Ensure Only');
title('Mean Distances');
grid on;

%% Correlate each neuron with each animals distance traveled
% For ethanol
numNeurons = size(animalNumSpikesEt, 1);
correlDistEt = [];

for iNeuron = 1:numNeurons
    neuronData = animalNumSpikesEt(iNeuron, 1:210);
    animalID = animalNumSpikesEt(iNeuron, 211);
    
    % Find the corresponding distances data for the current animal
    matchingRow = find(distEt(:, 211) == animalID, 1);  % Only take the first match
    
    % Only proceed if there is a match - rat 16 doesn't have any tracking
    % recorded so needed to include this so it skips 16
    if ~isempty(matchingRow)
        distData = distEt(matchingRow, 1:210);

        % Calculate correlation coefficient and p-value
        [corrCoeff, pValue] = corrcoef(neuronData, distData);
            
        % Append correls
        correlDistEt = [correlDistEt; corrCoeff(1, 2), pValue(1, 2)];
    end
end

% Same but for ensure
numNeurons = size(animalNumSpikesEn, 1);
numAnimals = size(distEn, 1);

correlDistEn = []; 
for iNeuron = 1:numNeurons
    neuronData = animalNumSpikesEn(iNeuron, 1:210);
    animalID = animalNumSpikesEn(iNeuron, 211);
    
    % Find the corresponding distances data for the current animal
    matchingRow = find(distEn(:, 211) == animalID, 1);  % Only take the first match
    
    % Only proceed if there is a match
    if ~isempty(matchingRow)
        distData = distEn(matchingRow, 1:210);

        % Calculate correlation coefficient and p-value
        [corrCoeff, pValue] = corrcoef(neuronData, distData);
            
        % Append correls
        correlDistEn = [correlDistEn; corrCoeff(1, 2), pValue(1, 2)];
    end
end
%% Find signficant correlations
%Find and index signficant positive and negative correlations for Ensure
alpha = 0.05;
sigPDistEn = correlDistEn(:, 2) < alpha;
sigR2DistEn = correlDistEn(sigPEn, 1);

posIndDistEn = sigR2DistEn > 0;
negIndDistEn = sigR2DistEn < 0;
totalDistInd = length(correlDistEn);
%Determine percents from the total
percentPosDistInd = sum(posIndDistEn) / totalDistInd * 100;
percentNegDistInd = sum(negIndDistEn) / totalDistInd * 100;
percentNonRespDist = 100 - (percentPosDistInd + percentNegDistInd);
percentPosNegEnDist = [percentPosDistInd, percentNegDistInd, percentNonRespDist];
%Graph in a pie chart
figure;
subplot(1,2,1);
pie(percentPosNegEnDist);
legend({'Positive', 'Negative', 'No Correl'}, 'Location', 'Best');
title('Correlations for Ensure Only')

%Same as above but with Ethanol
sigPDistEt = correlDistEt(:, 2) < alpha;
sigR2DistEt = correlDistEt(sigPEn, 1);

posIndDistEt = sigR2DistEt > 0;
negIndDistEt = sigR2DistEt < 0;
totalDistInd = length(correlDistEt);
%Determine percents from the total
percentPosDistInd = sum(posIndDistEt) / totalDistInd * 100;
percentNegDistInd = sum(negIndDistEt) / totalDistInd * 100;
percentNonRespDist = 100 - (percentPosDistInd + percentNegDistInd);
percentPosNegEtDist = [percentPosDistInd, percentNegDistInd, percentNonRespDist];
%Graph in a pie chart
subplot(1,2,2);
pie(percentPosNegEtDist);
legend({'Positive', 'Negative', 'No Correl'}, 'Location', 'Best');
title('Correlations for Ensure + EtOH')
%% Function for fill_between
function fill_between(x, y1, y2, color, alpha)  
    x = [x, fliplr(x)];
    y = [y1, fliplr(y2)];
    fill(x, y, color, 'FaceAlpha', alpha, 'EdgeAlpha', 0);
end