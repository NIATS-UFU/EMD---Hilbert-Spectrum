% Function name....: EstimateTFplots2
% Date.............: October 28, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This function estimates the Hilbert Spectrum, spectrogram and scalogram
%                    for a given signal (EMG). Contour plots are used for clarity.
%                    
% Parameters.......: 
%                    EMG  ........-> signal to be analysed
%                    fs   ........-> sampling frequency(Hz)
%                    NLevels .....-> number of contour plot levels
%                    tresholdHS...-> threshold for the HS
%                    tresholdSpec.-> threshold for the spectrogram
%                    tresholdScal.-> threshold for the scalogram
%                    tresholdWVD..-> threshold for the WVD

% Return...........:
%                    This function does not return any variable

function [] = EstimateTFplots(EMG,fs,NLevels,tresholdHS,tresholdSpec,tresholdScal,tresholdWVD)

t = [1:1:length(EMG)]/fs;

%Plotting the raw EMG signal 
figure;
subplot(5,1,1); plot(t,1e2*EMG);set(gca,'Xlim',[t(1) t(end)]);ylabel('\muV'); 
set(gca,'XtickLabel','');
drawnow; title('EMG signal','FontSize',9,'VerticalAlignment','middle');

%Estimating intrinsic mode functions
[IMFs,residue]= sig_to_imf(EMG',1e-5,2);

%Estimating Hilbert spectrum
n_bins = 400;
LFreq = 1; %1 minimum frequency in Hz
UFreq = fs/2; %maximum frequency in Hz
VisualizationMode = 1; % 2-D countour plot
subplot(5,1,2);
[m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(IMFs,t,fs,n_bins,[LFreq UFreq],VisualizationMode);


%Generating auxiliary time and frequency vectors
f = 1:1:n_bins;
f = convScale(1,n_bins,f,minFreq,maxFreq);

time = 1:1:n_bins;
time = convScale(1,n_bins,time,min(t),max(t));

%Performing energy normalization (between 0 and 1)
[m_a_p]= convScale(min(min(m_a_p)),max(max(m_a_p)),m_a_p,0,1);

%Estimating contour levels
CH = contourc(time,f,m_a_p,NLevels); %Estimating contour levels
[contourLevelsH] = getContourLevels(CH); %Getting contour levels

%Selecting contour levels above a pre-defined threshold
selContours = find(contourLevelsH>=tresholdHS);

subplot(5,1,2); %Selecting axis
CH = contour(time,f,m_a_p,contourLevelsH(selContours)); %Estimating contour plots 

colormap hot;
set(gca, 'color', 'k');
set(gca,'Ylim',[min(get(gca,'YLim')) 1000]);
title('');
set(gca,'XtickLabel','');
xlabel(''); ylabel('Hz'); title('Hilbert Spectrum','FontSize',9,'VerticalAlignment','middle');
set(gca,'YLim',[0 max(get(gca,'YLim'))]);

drawnow;

%%%%%%%%%%% Spectrogram 

[B_F,f_F,t_F] = specgram(EMG,256,fs,32); colormap hot;set(gca,'Color','k');

B_F =abs(B_F);% 20*log10(abs(B_F));%Getting spectrogram magnitude
%Performing energy normalization (between 0 and 1)
[B_F]= convScale(min(min(B_F)),max(max(B_F)),B_F,0,1);

%Estimating contour levels
CF = contourc(t_F,f_F,B_F,NLevels); %Estimating contour levels
[contourLevelsF] = getContourLevels(CF); %Getting contour levels

%Selecting contour levels above a pre-defined threshold
selContours = find(contourLevelsF>=tresholdSpec);

subplot(5,1,3); 
CF = contour(t_F,f_F,B_F,contourLevelsF(selContours));
colormap hot;set(gca,'Color','k') ;
ylabel('Hz'); xlabel('time (s)');

set(gca,'Ylim',[min(get(gca,'YLim')) 1000]);

set(gca,'XtickLabel','');
xlabel('');title('Spectrogram','FontSize',9,'VerticalAlignment','middle');
set(gca,'YLim',[0 max(get(gca,'YLim'))]);

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Scalogram

%Decomposing EMG signal
scales = [1:1:400];
c = cwt(EMG,scales,'db4');

%Performing energy normalization (between 0 and 1)
[Scal]= convScale(min(min(c)),max(max(c)),c,0,1);

%Converting scales into frequencies
Fwav = scal2frq(scales,'db4',1/fs);

%Estimating contour levels
CScal = contourc((1:length(EMG))/fs,Fwav,Scal,NLevels); %Estimating contour levels
[contourLevelsScal] = getContourLevels(CScal); %Getting contour levels

%Selecting contour levels above a pre-defined threshold
selContours = find(contourLevelsScal>=tresholdScal);

%Plotting
subplot(5,1,4); 
[Cwav] = contour((1:length(EMG))/fs,Fwav,Scal,contourLevelsScal(selContours));
colormap hot;set(gca,'Color','k');
set(gca,'XtickLabel','');
xlabel('');ylabel('Hz'); set(gca,'Ylim',[min(get(gca,'YLim')) 1000]);
title('Scalogram','FontSize',9,'VerticalAlignment','middle');
set(gca,'YLim',[0 max(get(gca,'YLim'))]);

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Wigner-Ville Distribution
%tfrwv is a member of the Time-Frequency Toolbox for Use with MatLab (Auger et al., CNRS France, 1996)
% 256 = number of frequency bins
% 1 = flag. Non zero if the status of the computation of the WVD is to be presented
%Note that an analystical signal instead of the real EMG signal is used for elimination
%or reduction of the interference of cross-terms
nbins = 256;
[TFR,T,F]=tfrwv(hilbert(EMG),[1:length(EMG)],nbins,1);
Freq = round(F*nbins/0.5); 
df = fs/nbins;%Frequency resolution (Hz)
Freq = Freq * df; %Frequency (Hz)

%Performing energy normalization (between 0 and 1)
[TFR]= convScale(min(min(TFR)),max(max(TFR)),TFR,0,1);

%Estimating contour levels
CTFR = contourc((1:length(EMG))/fs,Freq,TFR,NLevels); %Estimating contour levels
[contourLevelsCTFR] = getContourLevels(CTFR); %Getting contour levels

%Selecting contour levels above a pre-defined threshold
selContours = find(contourLevelsCTFR>=tresholdWVD);

%Displaying results
subplot(5,1,5);
contour((1:length(EMG))/fs,Freq,TFR,contourLevelsCTFR(selContours));

colormap hot;set(gca,'Color','k');
xlabel('time (s)'); ylabel('Hz'); set(gca,'Ylim',[min(get(gca,'YLim')) 1000]);
title('Wigner-Ville distribution','FontSize',9,'VerticalAlignment','middle');
set(gca,'YLim',[0 max(get(gca,'YLim'))]);
xlabel('time (s)');
drawnow;


%%%%%%%% Some final adjustments for data presentation
subplot(5,1,5);%Making the axis of the 4th figure the current one
Hfig5 = gca;
PosFigure5 = get(gca,'Position'); %Getting colorbar positioning and dimensions of the 4th figure
subplot(5,1,2);%Making the axis of the 4th figure the current one
PosFigure2 = get(gca,'Position'); %Getting colorbar positioning and dimensions of the 2nd figure

subplot(5,1,5);%Making the axis of the 4th figure the current one
hColorBar = colorbar; %displaying colorbar 
PosColorBar = get(hColorBar,'Position'); %Getting colorbar positioning and dimensions
set(hColorBar,'Position',[PosFigure5(1)+1.01*PosFigure5(3) PosFigure5(2)...
                         0.5*PosColorBar(3) PosFigure2(2)+PosFigure2(4)-PosFigure5(2)]); %PosFigure2(4)
set(hColorBar,'FontSize',9);

set(Hfig5,'Position',PosFigure5);%Restauring dimensions of subplot4
text(PosFigure5(1)+1.23*PosFigure5(3),abs(PosFigure5(2)-3.4),'normalized energy','Rotation',-90,'HorizontalAlignment','left',...
     'VerticalAlignment','bottom','Units',' normalized');