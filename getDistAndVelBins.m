clear all
close all

baseDir = uigetdir; 
folders = dir(baseDir);
folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));

allAnimalsDistData = struct();

% Constant values
xMin = -4.2; xMax = -0.6; % max and min x V
yMin = -5; yMax = 0; % max and min y V
numBins = 210; % number of bins to match neural data

allBinnedDistance = [];
allBinnedVelocity = [];

for f = 1:length(folders)
    animalID = folders(f).name;
    animalID_clean = matlab.lang.makeValidName(animalID); 
    animalPath = fullfile(baseDir, animalID);
    load(fullfile(animalPath, "stream.mat"));
    
    LFPs = double(traces) .* [channels.gain_to_uV] + [channels.offset_to_uV]; % convert from arbitrary units to uV
    
    % One file didn't record ADC and AUX channels, so skip if they dont
    % exist
    if size(LFPs, 2) < 70
        warning('Skipping %s: Missing columns 68 or 70 in LFPs.', animalID);
        continue;
    end
    
    time = time';
    anyMazeXY = [(LFPs(:,68)),(LFPs(:,70))];% columns corresponding to AUX1 and AUX3 for anyMaze tracking
    anyMazeXY = [time, anyMazeXY];
    anyMazeXY = anyMazeXY(1:1000/2:end, :);% down sample from 1000hz to 2hz

    time = anyMazeXY(:, 1);
    xuV = anyMazeXY(:, 2); % AUX1 = x coord in uV
    yuV = anyMazeXY(:, 3); % AUX3 = y coord in uV

    % Convert from microvolts (uV) to volts (V)
    x = xuV * 1e-6;
    y = yuV * 1e-6;

    % Clean up weird values outside of actual box by replacing with last in-bounds value
    validX = x;
    validY = y;

    for i = 2:length(x)
        if x(i) < xMin || x(i) > xMax
            validX(i) = validX(i-1);
        end
        if y(i) < yMin || y(i) > yMax
            validY(i) = validY(i-1);
        end
    end

    % Convert to inches (9.5 x 16.5" box)
    xInches = (validX - xMin) / (xMax - xMin) * 16.5;
    yInches = (validY - yMin) / (yMax - yMin) * 9.5;

    % Convert to cm, cause science
    xMM = xInches * 2.54;
    yMM = yInches * 2.54;

    % Using 210 bins to match the FR neural data
    binEdges = linspace(time(1), time(end), numBins + 1);
    binnedDistance = zeros(numBins, 1);
    binnedVelocity = zeros(numBins, 1);

    % Loop through bins and calculate distance and velocity
    for i = 1:numBins
        inBin = time >= binEdges(i) & time < binEdges(i + 1);
        
        % Get the x, y, and time values for the bin
        binX = xMM(inBin);
        binY = yMM(inBin);
        binTime = time(inBin);
        
        % Calculate distance traveled in bin
        if length(binX) > 1
            dx = diff(binX);
            dy = diff(binY);
            distances = sqrt(dx.^2 + dy.^2);
            binnedDistance(i) = sum(distances);
            
            % Calculate velocity: total bin distance / total bin time
            totalTime = binTime(end) - binTime(1);
            if totalTime > 0
                binnedVelocity(i) = binnedDistance(i) / totalTime;
            else
                binnedVelocity(i) = NaN;
            end
        else
            % No movement if point doesnt change so set to 0
            binnedDistance(i) = 0;
            binnedVelocity(i) = NaN;
        end
    end
    % Save with session info attached vs. not attached
    allAnimalsDistData.(animalID_clean).binnedDistance = binnedDistance;
    allAnimalsDistData.(animalID_clean).binnedVelocity = binnedVelocity;
    allBinnedDistance = [allBinnedDistance, binnedDistance];
    allBinnedVelocity = [allBinnedVelocity, binnedVelocity];
end
% Attach animal IDs based on order of sessions
labels = {'16en', '17en', '17et', '18et', '19en', 'NaN', '19et', '37et', '37en', 'NaN', '38et'};
allBinnedDistance = [num2cell(allBinnedDistance);labels];
allBinnedVelocity = [num2cell(allBinnedVelocity);labels];
% Export as .mat file
save('allAnimalsDistData.mat', 'allAnimalsDistData', 'allBinnedDistance', 'allBinnedVelocity');
