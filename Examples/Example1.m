

addpath( [pwd '\Examples'], [pwd '\Examples\SampleData'])

% Load EMG data
filename = 'SampleData.xlsx'; 
num = xlsread(filename, 'EMG data');

t = num(:,2);
EMG = num(:,3)';
plot(t, EMG);

% Call the DSP toolbox (set the sampling frequency to 5025 Hz in the
% interface sampling frequency field)
guiDSP