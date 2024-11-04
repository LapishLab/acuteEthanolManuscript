%clear all
%close all

enCumul = xlsread('ephysensurenoncumul.xlsx');
etCumul = xlsread('ephysetohnoncumul.xlsx');

%% Interpolate ensure and ethanol intakes
itpEnIntake = [];
for row = 1:size(enCumul, 1)
    x = enCumul(row, :);
    rowValues = [];
    for i = 1:size(x, 2)
        rowValues = [rowValues repmat(x(i)/10, 1, 10)];
    end
    itpEnIntake = [itpEnIntake; rowValues];
end

itpEtIntake = [];
for row = 1:size(etCumul, 1)
    x = etCumul(row, :);
    rowValues = [];
    for i = 1:size(x, 2)
        rowValues = [rowValues repmat(x(i)/10, 1, 10)];
    end
    itpEtIntake = [itpEtIntake; rowValues];
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