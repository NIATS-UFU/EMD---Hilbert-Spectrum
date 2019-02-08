% # EMD---Hilbert-Spectrum
% Matlab toolbox with several functions for time-frequency analysis
% 
% Contact: Prof. Adriano Andrade (adriano@ufu.br)
% 
% How to cite this toolbox:
% 
% [1] ANDRADE, A.; KYBERD, P.; NASUTO, S. The application of the Hilbert spectrum to the analysis of electromyographic signals. Inf. Sci. (Ny). 2008, 178, 2176–2193.
% https://doi.org/10.1016/j.ins.2007.12.013
% 
% [2] Andrade, A.O. (2005). Decomposition and analysis of electromyographic signals (PhD thesis). University of Reading, Reading, United Kingdom. 
% https://doi.org/10.13140/RG.2.2.33350.06727
% 
% https://ethos.bl.uk/OrderDetails.do?uin=uk.bl.ethos.422770

addpath( [pwd '\Examples'], [pwd '\Examples\SampleData'])

% Load EMG data
filename = 'SampleData.xlsx'; 
num = xlsread(filename, 'EMG data');

t = num(:,2); % time (s)
y = num(:,3)'; % EMG signal
plot(t, y);

fs = 5025; %sampling frequency


figure;
subplot(3,1,1);
plot(t,(y/10000)*1e6); xlabel('time (s)'); ylabel('\muV');

subplot(3,1,2);
[IMFs,residue]= sig_to_imf(y,1e-5,2); %Estimating intrinsic mode functions
LFreq=0;
UFreq=fs/2;
UFreq=500;
n_bins= 400;

[m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(IMFs,t,fs,n_bins,[LFreq UFreq],0); %Estimating Hilbert Spectrum

%Generating auxiliary time and frequency vectors
f = 1:1:n_bins;
f = convScale(1,n_bins,f,minFreq,maxFreq);
%Performing energy normalization (between 0 and 1)
[m_a_p]= convScale(min(min(m_a_p)),max(max(m_a_p)),m_a_p,0,1);


T=(0:1:n_bins-1)*hs_dt;
B=m_a_p;
F = f;
imagesc(T,f,m_a_p); axis xy; colormap('hot');drawnow;

[imnf_out2] = IMNF(y,B,F); %Estimating the Instantaneous mean frequency based on the Hilbert Spectrum

hold on; plot(T,imnf_out2,'g');ylabel('HS (Hz)'); xlabel('time (s)');


%Box plot of showing the distribution of the IMNF (the green signal on the
%joint time-frequency distribution
figure; boxplot(imnf_out2,'orientation','vertical');ylim([0 500]);ylabel('Instantaneous mean frequency (Hz)');set(gca,'XTickLabel','');
