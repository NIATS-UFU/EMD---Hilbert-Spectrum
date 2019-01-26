
load('CaioL5T1.txt');

EMG = CaioL5T1(:,2);

%Hilbert Spectrum
%generation of a chirp signal
fs = 10040; %sampling frequency (Hz)
t = [0:1:fs-1]*1/fs;

iO=1;
iF=10040;

for qq=1:30,

    y = EMG(iO:iF)'; %one-second-signal window to be analyzed

   %Estimating intrinsic mode functions
    figure;
    [IMFs,residue]= sig_to_imf(y,1e-5,2);
    LFreq=0;
    UFreq=fs/2;
    UFreq=500;
    n_bins= 400;
    [m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(IMFs,t,fs,n_bins,[LFreq UFreq],0);

    %Generating auxiliary time and frequency vectors
    f = 1:1:n_bins;
    f = convScale(1,n_bins,f,minFreq,maxFreq);
    %Performing energy normalization (between 0 and 1)
    [m_a_p]= convScale(min(min(m_a_p)),max(max(m_a_p)),m_a_p,0,1);


    T=(0:1:n_bins-1)*hs_dt;
    B=m_a_p;
    F = f;
    imagesc(T,f,m_a_p); axis xy; colormap('hot');drawnow;

    [imnf_out2] = IMNF(y,B,F);

    hold on; plot(T,imnf_out2,'w');
    %figure(1);subplot(1,5,2); boxplot(imnf_out2);ylim([0 500]);xlabel('Hilbert Spectrum');ylabel('Instantaneous mean frequency (Hz)');set(gca,'XTickLabel','');
    figure;
    boxplot(imnf_out2);
    resMed(qq) = med;%the variable med is exported to the Matlab worspace through the function boxplot (which I have modified for this purpose)
    close all;
    
    %updating indices
    iO=iF+1;
    iF=iO+fs-1;
end;