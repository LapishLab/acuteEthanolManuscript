%clear all
%close all

rnnPred = csvread('averagePredictionsEtOH.csv');
rnnPred = rnnPred';

%% Interpolate ensure and ethanol intakes
idxRnn = linspace(1, size(rnnPred, 1), size(rnnPred, 1));
idxitp = linspace(1, size(rnnPred, 1), 180);
itpRnnPred = zeros(180, size(rnnPred, 2));
for col = 1:size(rnnPred, 2)
    itpRnnPred(:, col) = interp1(idxRnn, rnnPred(:, col), idxitp, 'linear');
end
%% Plots to look at the interpolation vs. the raw intakes
%plot(itpEnIntake','b')
%legend('Interpolated intakes (180 bins)')
%figure;
%plot(intake,'r')
%legend('Original Intakes (18 bins)')

%% Extra stuff not needed for etohCorrel

%save('itpEsurIntake.mat','itpEnIntake')

%itpEnIntake(:, col) = interp1(idxintake, intake(:, col), idxitp, 'linear');