% Function name....: plotHS1
% Date.............: November 19, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    plotHS1 generates the Hilbert Spectrum of a signal.
% Parameters.......: 
%                    i_m_fs....-> intrinsic mode functions (matrix or
%                                 vector if only 1 component is available)
%                    t.........-> time vector
%                    fs........-> sampling frequency (Hz)
%                    bins .....-> number of bins for frequency axis 
%                    Freq .....-> vector containing the maximum and minimum frequencies. 
%                                 The first element is the min freq and second one is the max.
%                                 Case Freq = [] autoscale will be used
%                    type ....-> type of plot: 0 = pcolor, 1 = contour plot
%                                (2D) and 2 = contour plot (3D) 
% Return...........:
%                    m_a_p...............-> matrix with energy distribution
%                                           in dB
%                    minFreq.............-> minimum of Freq or minimum frequency
%                                           used in autoscale mode
%                    maxFreq.............-> maximum of Freq or maximum frequency
%                                           used in autoscale mode
%                    hs_dt...............-> time resolution of the HS
%                    time................-> time vector
%                    f...................-> frequency vector  


function [m_a_p,minFreq,maxFreq, hs_dt,time,f] = plotHS1(i_m_fs,t,fs,bins, Freq,type)

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
    
%Estimation of intantaneous atributes of each IMF
for k = 1:r,
    [amp_,phase_,freq_]= instantAtrib(i_m_fs(k,:),fs);
    amp(k,:) = amp_;
    phase(k,:) = phase_;
    freq(k,:) = freq_;
end;

%Estimating scale
if(isempty(Freq)==1),
    minFreq = min(min(freq));    
    maxFreq = max(max(freq));
else
    minFreq = Freq(1);
    maxFreq = Freq(2);
end;%if

%generating time vector (this is a scaled version based on the number of bins)
t_s = round(convScale(min(t),max(t),t,1,bins)); %color map column   

%Scaling instantaneous frequency (based on the number of bins)
ii=[];
for k = 1:r,
    if(min(freq(k,:))<minFreq | max(freq(k,:))>maxFreq),
        ii = find(freq(k,:)>maxFreq);
        freq(k,ii) = minFreq;
        w_s(k,:) = round(convScale(minFreq,maxFreq,freq(k,:),1,bins)); %color map row
    else
        w_s(k,:) = round(convScale(minFreq,maxFreq,freq(k,:),1,bins)); %color map row
    end;
end;
    
%creating color map (Converting scale to dB)
m_a_p = 10*log10(min(min(amp)))*ones(bins);
    
for k = 1:r,
    for i=1:c,
        m_a_p(w_s(k,i),t_s(i)) = 10*log10(amp(k,i));
    end;
end;
   
%generation of frequency axis
f = 1:1:bins;
f = convScale(1,bins,f,minFreq,maxFreq);

%generating time axis
time = 1:1:bins;
time = convScale(1,bins,time,min(t),max(t));
hs_dt = time(2)-time(1);%time resolution of the HS
    
    
%%%%%%%%%%%%%%%%%%Plotting
switch (type)
    case 1,
        contour(time,f,m_a_p);
        title('\it{Hilbert spectrum - contourplot}','FontSize',10)
    case 2,
        contour3(time,f,m_a_p);
        title('\it{Hilbert spectrum} - countourplot/slice','FontSize',10)
    otherwise,
        pcolor(time,f,m_a_p);
        
        %imagesc(time,f,m_a_p);axis xy;
        shading('interp');
        title('\it{Hilbert spectrum}','FontSize',10)
end;

colormap hot;
colorbar('horiz');
ylabel('frequency (Hz)','FontSize',10);
xlabel('time (s)','FontSize',10);