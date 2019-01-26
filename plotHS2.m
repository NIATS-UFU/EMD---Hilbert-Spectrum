% Function name....: plotHS2
% Date.............: November 19, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    plotHS2 generates the Hilbert Spectrum of a signal.
% Parameters.......: 
%                    i_m_fs....-> intrinsic mode functions (matrix or
%                                 vector if only 1 component is available)
%                    fs........-> sampling frequency (Hz)
%                    bins .....-> number of bins for frequency axis 
%                    Freq .....-> vector containing the maximum and minimum frequencies. 
%                                 The first element is the min freq and second one is the max.
%                                 Case Freq = [] autoscale will be used
% Return...........:
%                    m_a_p...............-> matrix with energy distribution
%                    minFreq.............-> minimum of Freq or minimum frequency
%                                           used in autoscale mode
%                    maxFreq.............-> maximum of Freq or maximum frequency
%                                           used in autoscale mode
%                    normMin.............-> the matrix m_a_p is normalized
%                                           between 0 and 1. For recovering
%                                           the original scale in dB
%                                           normMin should be used as the
%                                           minimum limit of the original
%                                           scale
%                    normMax.............-> the matrix m_a_p is normalized
%                                           between 0 and 1. For recovering
%                                           the original scale in dB
%                                           normMax should be used as the
%                                           maximum limit of the original
%                                           scale
% Example..........:
%                   1) linear chirp signal
%                      t = 0:0.001:2;  % 2 secs @ 1kHz sample rate
%                      y = chirp(t,0,1,150); % Start @ DC, cross 150Hz at t=1 sec
%                      plotHS2(y,1000,400,[0 500]);
%                   2) quadratic chirp signal 
%                      t = -2:0.001:2; % ±2 secs @ 1kHz sample rate
%                      y = chirp(t,100,1,200,'quadratic'); 
%                      plotHS2(y,1000,400,[0 500]);


function [m_a_p,minFreq,maxFreq,normMin,normMax] = plotHS2(i_m_fs,fs,bins,Freq)

%Getting dimensions of the matrix (or vector) i_m_fs
[r,c]=size(i_m_fs); % r= number of imfs

%We do want to work with information stored in rows
if r>c,
    i_m_fs = i_m_fs';
    [r,c]=size(i_m_fs); 
end;
    
%initialization of matrices
amp = zeros(r,c);  %intantaneous amplitude
phase = zeros(r,c); %instantaneous phase
freq = zeros(r,c); %instantaneous frequency
w_s = zeros(r,c); %scaled version of freq
    
%creating time vector
t = [0:1:c-1]/fs;
        
%Estimation of intantaneous atributes of each IMF
for k = 1:r,
    [amp_,phase_,freq_]= instantAtrib(i_m_fs(k,:),fs);
    amp(k,:) = amp_;
    phase(k,:) = phase_;
    freq(k,:) = freq_;
end;

%Estimating frequency scale range (it may be user defined)
if(isempty(Freq)==1),
    minFreq = min(min(freq));    
    maxFreq = max(max(freq));
else
    minFreq = Freq(1);
    maxFreq = Freq(2);
end;%if
    
%Scaling instantaneous frequency (based on the number of bins)
ii=[];
for k = 1:r,
    if(min(freq(k,:))<minFreq | max(freq(k,:))>maxFreq),
        ii = find(freq(k,:)>maxFreq);
        freq(k,ii) = minFreq;
        w_s(k,:) = round(convScale(minFreq,maxFreq,freq(k,:),1,bins)); %color map row
    else
        w_s(k,:) = round(convScale(minFreq,maxFreq,freq(k,:),1,bins)); %color map row
    end; %if
end;%for
    
%creating color map (Converting scale to dB)
m_a_p = 20*log10(min(min(amp)))*ones(bins,size(amp,2));

for k = 1:r,
    for i=1:size(amp,2),
        m_a_p(w_s(k,i),i) = 20*log10(amp(k,i));
    end;
end;

%generation of frequency axis
f = 1:1:bins;
f = convScale(1,bins,f,minFreq,maxFreq);

%normalization (between 0 and 1)
normMin = min(min(m_a_p));
normMax = max(max(m_a_p));
m_a_p = convScale(normMin,normMax,m_a_p,0,1);

%%%%%%%%%%%%%%%%%%Plotting
imagesc(t,f,m_a_p);
colormap 'hot'; 
axis xy;
ylabel('frequency (Hz)','FontSize',10);
xlabel('time (s)','FontSize',10);