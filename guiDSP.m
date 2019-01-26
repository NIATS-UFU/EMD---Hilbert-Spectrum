function varargout = guiDSP(varargin)
% GUIDSP Application M-file for guiDSP.fig
%    FIG = GUIDSP launch guiDSP GUI.
%    GUIDSP('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 09-Dec-2014 21:26:39

if nargin == 0  % LAUNCH GUI
	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.




% --------------------------------------------------------------------

function varargout = edit2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit5_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkOffset_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;
global offset;

N = length(s_i_g_n_a_l);

sampFreq = str2double(get(handles.editSampFreq,'String'));
%axes(handles.axesSignal);
t=0:1/sampFreq:(N-1)/sampFreq;

state = get(handles.checkOffset,'Value');

if state==1,
    offset = mean(s_i_g_n_a_l);
    s_i_g_n_a_l = s_i_g_n_a_l-offset;
     plot(t,s_i_g_n_a_l);
else
    s_i_g_n_a_l = s_i_g_n_a_l+offset;
    plot(t,s_i_g_n_a_l);
end; %if

% --------------------------------------------------------------------
function varargout = checkOffset_CreateFcn(h, eventdata, handles, varargin)

% global s_i_g_n_a_l;
% 
% N = length(s_i_g_n_a_l);
% 
% sampFreq = str2double(get(handles.editSampFreq,'String'));
% axes(handles.axesSignal);
% t=0:1/sampFreq:(N-1)/sampFreq;
% plot(t,s_i_g_n_a_l);




% --------------------------------------------------------------------
function varargout = axesSignal_CreateFcn(h, eventdata, handles, varargin)

% global s_i_g_n_a_l;
% global offset;
% 
% N = length(s_i_g_n_a_l);
% 
% %sampFreq = str2double(get(handles.editSampFreq,'String'));
% %axes(handles.axesSignal);
% sampFreq = 4000;
% t=0:1/sampFreq:(N-1)/sampFreq;
% 
% plot(t,s_i_g_n_a_l);
% 
% % if state==1,
% %     offset = mean(s_i_g_n_a_l);
% %     s_i_g_n_a_l = s_i_g_n_a_l-offset;
% %     plot(t,s_i_g_n_a_l);
% % else
% %     s_i_g_n_a_l = s_i_g_n_a_l+offset;
% %     plot(t,s_i_g_n_a_l);
% % end; %if


% --------------------------------------------------------------------
function varargout = figure1_CreateFcn(h, eventdata, handles, varargin)

vars = evalin('base','who');
disp(handles);
% set(handles.listbox_var,'String','vars');


% --------------------------------------------------------------------
function varargout = menu_dsp_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = menu_extrema_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
  
    %calculating signal extrema (valleys and peaks).
    [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(s_i_g_n_a_l);
    
    %exporting instantaneou attributes to workspace
    assignin('base','LocalMinAmplitude',LocalMinV);
    assignin('base','LocalMaxAmplitude',LocalMaxV);
    assignin('base','LocalMinIndex',LocalMinI);
    assignin('base','LocalMaxIndex',LocalMaxI);

    figure(1);
    plot(s_i_g_n_a_l);
    hold on; plot(LocalMaxI,LocalMaxV,'r+');
    hold on; plot(LocalMinI,LocalMinV,'y*');
    legend('signal','peaks','valleys');

% --------------------------------------------------------------------
function varargout = menu_hist_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
    
    Nbits = str2double(get(handles.Nbits,'String'));

    %Estimating the histogram
    [xout,n]= histogram(s_i_g_n_a_l,Nbits);
 
    figure(1);
    title('Histogram');
    plot(xout,n);
    ylabel('Number of ocurrences');
    xlabel('Amplitude (V)');

% --------------------------------------------------------------------
function varargout = Untitled_4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = menu_env_Callback(h, eventdata, handles, varargin)

  global s_i_g_n_a_l;

    %Getting the interpolation method
    index_selected = get(handles.edit_interpMethod,'Value');

    %Estimating the signal envelope
    [AveEnvelope,yu,yl,n_extrema]= meanEnv(s_i_g_n_a_l,index_selected-1);
    
    %exporting envelopes to workspace
    assignin('base','UpperEnvelope',yu);
    assignin('base','LowerEnvelope',yl);
    assignin('base','AveEnvelope',AveEnvelope);

    %Creating figure
    handle_figure = 333;
    figure(handle_figure);
    set(handle_figure,'Name','Envelope');
    set(handle_figure,'NumberTitle','off');

    %Plotting
    plot(s_i_g_n_a_l);
    hold on; plot(yu,'r--');
    hold on; plot(yl,'y--');
    hold on; plot(AveEnvelope,'k','LineWidth',2);
    
    ylabel('Amplitude (V)');
    xlabel('Samples');
    legend('signal','upper envelope','lower envelope', 'mean envelope');

% --------------------------------------------------------------------
function varargout = menu_burst_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
    
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    maxADRange = str2double(get(handles.maxADRange,'String'));
    Nbits = str2double(get(handles.Nbits,'String'));
    
    [onset,tr]= detectBurst(s_i_g_n_a_l,maxADRange,Nbits,sampFreq);
    figure(1);
    plot(onset);
    hold on; plot(s_i_g_n_a_l,'r');
    hold on; plot([0 length(s_i_g_n_a_l)],[tr tr],'g','LineWidth',2);
    hold on; plot([0 length(s_i_g_n_a_l)],[-tr -tr],'g','LineWidth',2); 
    ylabel('Amplitude (volts)');
    xlabel('Samples');
    
    %exporting onset to workspace
    assignin('base','Onset',tr);


% --------------------------------------------------------------------
function varargout = menu_InstAt_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;

   
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    
    [amp,phase,freq]= instantAtrib(s_i_g_n_a_l,sampFreq);
    
    %exporting instantaneou attributes to workspace
    assignin('base','amp',amp);
    assignin('base','phase',phase);
    assignin('base','freq',freq);
    
    %Plotting
    handle_figure = 222;
    figure(handle_figure);
    set(handle_figure,'Name','Intantaneous atributes');
    set(handle_figure,'NumberTitle','off');
    
    subplot(3,1,1);plot(amp);
    ylabel('amplitude (V)');

    subplot(3,1,2);plot(phase);
    ylabel('phase (rad)');
    
    subplot(3,1,3);plot(freq);
    ylabel('frequency (Hz)');
    xlabel('samples');


% --------------------------------------------------------------------
function varargout = menu_psd_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;

sampFreq = str2double(get(handles.editSampFreq,'String'));

% estimates the Power Spectral Density of 
% a discrete-time signal vector X using Welch's averaged, modified
% periodogram method.  
[Pxx,F] = PSD(s_i_g_n_a_l,256,sampFreq);
Pxx = 10*log10(Pxx); %dB

%Plotting
handle_figure = 4;
figure(handle_figure);
set(handle_figure,'Name','PSD - based on Welch`s averaged, modified periodogram method');
set(handle_figure,'NumberTitle','off');
    
Pxx_norm = convScale(min(Pxx),max(Pxx),Pxx,0,1); %normalizing the spectrum [0 1]
plot(F,Pxx_norm);
xlabel('frequency (Hz)');
ylabel('Normalized amplitude');
title('Power Spectral Density estimate');

% --------------------------------------------------------------------
function varargout = menu_emd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = menu_sp_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
    global i_m_fs;
    global residue;
    
    i_m_fs = [];
    
    %Getting mse (defined by the user)
    mse = str2double(get(handles.edit_mse,'String'));
    
   
    %Getting the interpolation method
    index_selected = get(handles.edit_interpMethod,'Value');

    %Sifting process
    [i_m_fs,residue]= sig_to_imf(s_i_g_n_a_l,mse,index_selected-1);
    [r,c] = size(i_m_fs);
    
    %sampling frequency (defined by the user)
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    %Estimating time
    n=length(s_i_g_n_a_l);
    t=0:1/sampFreq:(n-1)/sampFreq;
    
    %Plotting
    handle_figure = 5;
    figure(handle_figure);
    set(handle_figure,'Name','Intrinsic mode functions');
    set(handle_figure,'NumberTitle','off');

    for i=1:r,
        subplot(r,1,i);plot(t,i_m_fs(i,:));
        if(i~=r),
            set(gca,'XTickLabel',{''});
        end;
    end%for
    xlabel('Time (s)');

%     %exporting imfs to workspace
    for i=1:r,
         assignin('base',['imf',num2str(i)],i_m_fs(i,:));
    end%for
    
    
% --------------------------------------------------------------------
function varargout = menu_hs_Callback(h, eventdata, handles, varargin)

    
  
    
    
% --------------------------------------------------------------------
function varargout = menu_mhs_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
    global m_a_p;
    global hs_dt;
    
%    Verifying if there exists a valid map to generate the Hilbert spectrum
    if isempty(m_a_p)==1,
        msgbox('There is Hilbert spectrum available.','guiDSP','error');
        return;
    end
    
    %Getting the lower frequency
    LFreq = str2double(get(handles.edit_minFreq,'String'));
    %Getting the upper frequency
    UFreq =str2double(get(handles.edit_maxFreq,'String'));  
    %Getting the sampling frequency
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    
    %Estimating the Marginal Hilbert spectrum
    [mhs,f] = MHSpec(m_a_p,LFreq,UFreq,hs_dt);
   
    %Plotting...
    handle_figure = 6;
    figure(handle_figure);
    set(handle_figure,'Name','Marginal Hilbert spectrum');
    set(handle_figure,'NumberTitle','off');
    mhs_norm = convScale(min(mhs),max(mhs),mhs,0,1); %normalizing the spectrum [0 1]
    plot(f,mhs_norm);
    xlabel('frequency (Hz)');
    ylabel('Normalized amplitude');

% --------------------------------------------------------------------
function varargout = menu_ds_Callback(h, eventdata, handles, varargin)

   global s_i_g_n_a_l;
   global m_a_p;
   global hs_dt; 
   
   if isempty(m_a_p)==1,
        msgbox('There is Hilbert spectrum available.','guiDSP','error');
        return;
   end
   
    %Getting the lower frequency
    LFreq = str2double(get(handles.edit_minFreq,'String'));
    %Getting the upper frequency
    UFreq =str2double(get(handles.edit_maxFreq,'String'));  
    %Getting samplinf frequency
    sampFreq = str2double(get(handles.editSampFreq,'String'));
   
    T = length(s_i_g_n_a_l)*1/sampFreq;

    [meanMHS,d_s,f] = DS(m_a_p,LFreq,UFreq,hs_dt,T);
    
    %Plotting...
    handle_figure = 60;
    figure(handle_figure);
    set(handle_figure,'Name','Degree of stationarity');
    set(handle_figure,'NumberTitle','off');
    figure(handle_figure);
    plot(f,d_s);
    xlabel('frequency');
    ylabel('DS');

% --------------------------------------------------------------------
function varargout = R_M_S_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;

prompt = {'Enter the window size'};
dlg_title = 'RMS estimation';
num_lines= 1;
def = {'20'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);
     
if(isempty(msgbox_output)==1),
    return;
else
    %Defining window size
    N = str2double(msgbox_output);
    %Root mean square estimation
    [rms,index] = rmsEst(s_i_g_n_a_l,N);

    handle_figure = 1111;
    figure(handle_figure);
    set(handle_figure,'Name','RMS');
    set(handle_figure,'NumberTitle','off');

    figure(handle_figure);

    plot(index,rms);
    hold on; plot(s_i_g_n_a_l,'r');

    %exporting filtered signal to workspace
    assignin('base','r_m_s',rms);
    assignin('base','index',index);

end; %if

% --------------------------------------------------------------------
function varargout = menu_spectrogram_Callback(h, eventdata, handles, varargin)

   global s_i_g_n_a_l;

   %Getting sampling frequency
   sampFreq = str2double(get(handles.editSampFreq,'String'));

   %Plotting...
   handle_figure = 1000;
   figure(handle_figure);
   set(handle_figure,'Name','Spectrogram');
   set(handle_figure,'NumberTitle','off');

   specgram(s_i_g_n_a_l,256,sampFreq,32);

   colormap(hot(256));
   colorbar;
   ylabel('frequency (Hz)');
   xlabel('time (s)');


% --------------------------------------------------------------------
function varargout = pushbutton7_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = editSampFreq_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = edit6_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit7_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = listbox2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = listbox_var_CreateFcn(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = listbox_var_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;
global offset;

%showing the command window variables in the list box
vars = evalin('base','who');
set(handles.listbox_var,'String',vars);

%selecting the variable chosen by the user
index_selected = get(handles.listbox_var,'Value');
s = vars{index_selected};
s_i_g_n_a_l = evalin('base',s);
axes(handles.axesSignal);%Setting the current axes
sampFreq = str2double(get(handles.editSampFreq,'String'));

%Estimating time
n=length(s_i_g_n_a_l);
t=0:1/sampFreq:(n-1)/sampFreq;
plot(t,s_i_g_n_a_l);
xlabel('Time (s)');

%updating offset
offset = 0;
set(handles.checkOffset,'Value',0);
% --------------------------------------------------------------------
function varargout = edit_mse_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit_tag_Callback(h, eventdata, handles, varargin)

global n_bins;
n_bins = str2double(get(handles.editSampFreq,'String'));


% --------------------------------------------------------------------
function varargout = edit_tag_CreateFcn(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = minADRange_Callback(h, eventdata, handles, varargin)








% --------------------------------------------------------------------
function varargout = edit10_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = edit_interpMethod_Callback(h, eventdata, handles, varargin)





% --------------------------------------------------------------------
function varargout = IMF_anlysis_Callback(h, eventdata, handles, varargin)

    global s_i_g_n_a_l;
    global i_m_fs;

    Pxx=[];
    F = [];
    str = [];
    
    %Getting sampling frequency
    sampFreq = str2double(get(handles.editSampFreq,'String'));

    % estimates the Power Spectral Density of 
    % a discrete-time signal vector X using Welch's averaged, modified
    % periodogram method.  
    WindowSize = 256;
    [Pxx,F] = PSD(s_i_g_n_a_l,WindowSize,sampFreq); %input signal
    MaxPxx=max(Pxx);
    str  = ('signal');

%     %Estimating PSD function of IMFs
    numberIMFs = min(size(i_m_fs));
% 
     for i=1:numberIMFs,
           [Pxx(:,i+1),F(:,i+1)] = PSD(i_m_fs(i,:),WindowSize,sampFreq); %input signal
           auxStr= (['imf ', num2str(i)]);
           str = str2mat(str,auxStr);
     end;
% 
%     %Plotting
    handle_figure = 114;
    figure(handle_figure);
    set(handle_figure,'Name','IMF analysis');
    set(handle_figure,'NumberTitle','off');
    plot(F,Pxx/MaxPxx);
    legend(str);
    xlabel('Frequency (Hz)');
    ylabel('Normalized amplitude');
    title('Power spectrum');

% --------------------------------------------------------------------
function varargout = menu_completeness_Callback(h, eventdata, handles, varargin)

global i_m_fs;
global residue;
global s_i_g_n_a_l;

%Reconstructing data
aux=[i_m_fs;residue];
SumAux = sum(aux);
    
handle_figure = 200;
figure(handle_figure);
set(handle_figure,'Name','Completeness analysis');
set(handle_figure,'NumberTitle','off');
    
plot(s_i_g_n_a_l);
hold on; plot(SumAux,'r:');
hold on; plot(s_i_g_n_a_l-SumAux,'k-.');
legend('signal', 'reconstructed signal', 'error');
xlabel('Number of samples');
ylabel('Amplitude (V)');

% --------------------------------------------------------------------
function varargout = WLPD_filter_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;
global WLPDsignal;

%Getting sampling frequency
sampFreq = str2double(get(handles.editSampFreq,'String'));

%filtering
[WLPDsignal] = WLPD(s_i_g_n_a_l,sampFreq,40);

handle_figure = 211;
figure(handle_figure);
set(handle_figure,'Name','Weighted low pass differential filter');
set(handle_figure,'NumberTitle','off');
plot(WLPDsignal);

%exporting filtered signal to workspace
assignin('base','WLPDsignal',WLPDsignal);

% --------------------------------------------------------------------
function varargout = MUAP_Analysis_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = maxADRange_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = listbox_var_ButtonDownFcn(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = menu_Principal_peaks_Callback(h, eventdata, handles, varargin)

   global s_i_g_n_a_l;
    
   sampFreq = str2double(get(handles.editSampFreq,'String'));
   maxADRange = str2double(get(handles.maxADRange,'String'));
   Nbits = str2double(get(handles.Nbits,'String'));
    
   [Index,Amplitude]= mainPeaks(s_i_g_n_a_l,maxADRange,Nbits,sampFreq);

   figure(1);
   plot(s_i_g_n_a_l);
   hold on; plot(Index,Amplitude,'ko');




% --------------------------------------------------------------------
function varargout = Filter_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = MovAvgFilt_Callback(h, eventdata, handles, varargin)
    
    %calling filter order setup window
    MovAvgFilter;
    




% --------------------------------------------------------------------
function varargout = EMG_simulator_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = menu_EFAP_Callback(h, eventdata, handles, varargin)
 
%calling window
%GUI_AP;
[ep,z,t] = GenerateFAP;

figure(1);
plot(t,ep);
title('Muscle fibre action potential');
xlabel('time (ms)');
ylabel('amplitude (mV)');
%exporting filtered signal to workspace
assignin('base','ep',ep);
assignin('base','z',z);
assignin('base','t',t);
% --------------------------------------------------------------------
function varargout = menu_MUAP_Callback(h, eventdata, handles, varargin)

 %calling window
%GUI_MUAP;



% --------------------------------------------------------------------
function varargout = menu_MUAPT_Callback(h, eventdata, handles, varargin)

 %calling window
%GUI_MUAPT;


% --------------------------------------------------------------------
function varargout = menu_EMG_Callback(h, eventdata, handles, varargin)
%calling window
%GUI_EMG;

MU = evalin('base','MU');
Fs = evalin('base','Fs');

if(isempty(MU)==1 | isempty(Fs)==1  ),
     msgbox('MU/Fs is not defined.','guiDSP','error');
     return;
end;

prompt = {'SNR (dB)'};
dlg_title = 'EMG signal simulation';
num_lines= 1;
def = {'80'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    SNR = str2num(msgbox_output{1});
    [EMG]= SimulateCompEMGsignal(MU,Fs,SNR);
end;

figure(3333);
%Number of electrodes
timeVec = 1/Fs*[0:1:size(EMG,2)-1];
plot(timeVec,EMG');

%exporting filtered signal
assignin('base','EMG',EMG);

% --------------------------------------------------------------------
function LPButterworth_Callback(hObject, eventdata, handles)
% hObject    handle to LPButterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 %calling window
global s_i_g_n_a_l;

prompt = {'Cutoff frequency (Hz)', 'order'};
dlg_title = 'Low-pass Butterworth filter';
num_lines= 1;
def = {'20','4'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    cutoffFreq = str2num(msgbox_output{1});
    order = str2num(msgbox_output{2});
    %current sampling frequency 
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    Wn = cutoffFreq/(sampFreq/2);
    [b,a] = butter(order,Wn); %filter design
    figure(1); freqz(b,a,round(sampFreq/2),sampFreq); %plotting filter characteristics
    y = filtfilt(b,a,s_i_g_n_a_l); %zero-phase digital filtering (zero-phase distortion)
    %exporting filtered signal
    assignin('base','LowPassButFiltSig',y);
end;

%-------------------
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over editSampFreq.
function editSampFreq_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to editSampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function editSampFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSampFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --------------------------------------------------------------------
function ss_menu_nws_Callback(hObject, eventdata, handles)
% hObject    handle to ss_menu_nws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s_i_g_n_a_l;
global r;

if(isempty(s_i_g_n_a_l)==0),

    sampFreq = str2double(get(handles.editSampFreq,'String'));

    %selecting the variable chosen by the user
    index_selected = get(handles.listbox_var,'Value');
    s = get(handles.listbox_var,'String');
    s = s{index_selected};

    %Creating sptool (It is required to have the MatLab signal processing
    %toolbox for this option
    sptool('load','Signal',s_i_g_n_a_l,sampFreq,s);
end;

% --------------------------------------------------------------------
function ss_menu_peak_detector_Callback(hObject, eventdata, handles)
% hObject    handle to ss_menu_peak_detector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s_i_g_n_a_l;
% global r; %structure containing the selected markers (noise window, etc)

%guaranteeing that samples are in columns
  [r,c] = size(s_i_g_n_a_l);
  if (r>c),
      s_i_g_n_a_l = s_i_g_n_a_l';
  end;

%showing the command window variables in the list box
vars = evalin('base','who');
%selecting the variable chosen by the user
% index_selected = get(handles.listbox_var,'Value');
% s = vars{index_selected};
r = evalin('base','r');%structure containing the selected markers (noise window, etc)

if(isempty(r)==1),
     msgbox('A noise of window should be selected and its markers should be exported to a workspace variable called r .','guiDSP','error');
     return;
end;

if(isempty(s_i_g_n_a_l)==0),
    
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    [amp,phase,freq]= instantAtrib(s_i_g_n_a_l,sampFreq);% estimating instantaneous attributes of the signal
    
    noise = amp(round(r.x1*sampFreq):round(r.x2*sampFreq));

    th = 5*std(noise); %estimation of the noise based on the standard deviation of the chosen window
    
    %DetectingMUAPs
    [wpeak,wnd]= ss_peakDetector(amp,th,sampFreq);
    
    %Detected MUAPs to the workspace
    assignin('base','ss_wpeak',wpeak);
    assignin('base','ss_wnd',wnd);

    figure(1);
    plot(amp);
    hold on; stem(wpeak.pos, wpeak.amp,'r');
    hold on; plot(wnd.pos_i, wnd.amp_i,'k>');
    hold on; plot(wnd.pos_f, wnd.amp_f,'k<');
    hold on; plot([0 length(amp)],[th th],'g:');
 
    %Calling a dialog box for definition of the number of bins of the
    %histogram
    prompt = {'Enter the number of bins for the MUAP histogram:'};
    dlg_title = 'Histogram';
    num_lines= 1;
    msgbox_output = inputdlg(prompt,dlg_title,num_lines);
    
    if(isempty(msgbox_output)==1),
        return;
    else
      bins  = str2double(msgbox_output);
      figure(2);
      hist(wpeak.amp,bins);
      xlabel('peaks - amplitude (V)');
      ylabel('number of ocurrences');%frequency
    end;
end;%if


% --------------------------------------------------------------------
function HS_pcolour_Callback(hObject, eventdata, handles)
% hObject    handle to HS_pcolour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  global s_i_g_n_a_l;
    global i_m_fs;
    global m_a_p;
    global n_bins;
    global hs_dt;
    
    if isempty(i_m_fs)==1,
        msgbox('There is no intrinsic mode function available.','guiDSP','error');
        return;
    end

    %sampling frequency (defined by the user)
    sampFreq = str2double(get(handles.editSampFreq,'String'));
 
    %Getting the lower frequency
    LFreq = str2double(get(handles.edit_minFreq,'String'));
    %Getting the upper frequency
    UFreq =str2double(get(handles.edit_maxFreq,'String'));  
    
    %Getting the number of bins
    n_bins = str2double(get(handles.edit_bins,'String'));  
    
    %Estimating time
    n=length(s_i_g_n_a_l);
    t=0:1/sampFreq:(n-1)/sampFreq;
   
    handle_figure = 7;
    figure(handle_figure);
    set(handle_figure,'Name','Hilbert spectrum');
    set(handle_figure,'NumberTitle','off');

    subplot(2,1,1); plot(t,s_i_g_n_a_l);
    subplot(2,1,2); 
    [m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(i_m_fs,t,sampFreq,n_bins,[LFreq UFreq],0);

% --------------------------------------------------------------------
function HS_contourPlot_Callback(hObject, eventdata, handles)
% hObject    handle to HS_contourPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  global s_i_g_n_a_l;
    global i_m_fs;
    global m_a_p;
    global n_bins;
    global hs_dt;
    
    if isempty(i_m_fs)==1,
        msgbox('There is no intrinsic mode function available.','guiDSP','error');
        return;
    end

    %sampling frequency (defined by the user)
    sampFreq = str2double(get(handles.editSampFreq,'String'));
 
    %Getting the lower frequency
    LFreq = str2double(get(handles.edit_minFreq,'String'));
    %Getting the upper frequency
    UFreq =str2double(get(handles.edit_maxFreq,'String'));  
    
    %Getting the number of bins
    n_bins = str2double(get(handles.edit_bins,'String'));  
    
    %Estimating time
    n=length(s_i_g_n_a_l);
    t=0:1/sampFreq:(n-1)/sampFreq;
   
    handle_figure = 7;
    figure(handle_figure);
    set(handle_figure,'Name','Hilbert spectrum');
    set(handle_figure,'NumberTitle','off');

    subplot(2,1,1); plot(t,s_i_g_n_a_l);
    subplot(2,1,2); 
    
    [m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(i_m_fs,t,sampFreq,n_bins,[LFreq UFreq],1);

% --------------------------------------------------------------------
function HS_3DContourPlot_Callback(hObject, eventdata, handles)
% hObject    handle to HS_3DContourPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  global s_i_g_n_a_l;
  global i_m_fs;
  global m_a_p;
  global n_bins;
  global hs_dt;
    
    if isempty(i_m_fs)==1,
        msgbox('There is no intrinsic mode function available.','guiDSP','error');
        return;
    end

    %sampling frequency (defined by the user)
    sampFreq = str2double(get(handles.editSampFreq,'String'));
 
    %Getting the lower frequency
    LFreq = str2double(get(handles.edit_minFreq,'String'));
    %Getting the upper frequency
    UFreq =str2double(get(handles.edit_maxFreq,'String'));  
    
    %Getting the number of bins
    n_bins = str2double(get(handles.edit_bins,'String'));  
    
    %Estimating time
    n=length(s_i_g_n_a_l);
    t=0:1/sampFreq:(n-1)/sampFreq;
   
    handle_figure = 7;
    figure(handle_figure);
    set(handle_figure,'Name','Hilbert spectrum');
    set(handle_figure,'NumberTitle','off');


    [m_a_p,minFreq,maxFreq, hs_dt] = plotHS1(i_m_fs,t,sampFreq,n_bins,[LFreq UFreq],2);


% --------------------------------------------------------------------
function ResampleTS_Callback(hObject, eventdata, handles)
% hObject    handle to ResampleTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Calling a dialog box for definition of the new sampling frequency 

global s_i_g_n_a_l;

prompt = {'Enter the new sampling frequency (Hz):'};
dlg_title = 'Resampling time-series';
num_lines= 1;
msgbox_output = inputdlg(prompt,dlg_title,num_lines);
    
if(isempty(msgbox_output)==1),
    return;
else
 %current sampling frequency 
 sampFreq = str2double(get(handles.editSampFreq,'String'));

 %new sampling frequency defined by the user
 newSamplingFreq  = str2double(msgbox_output);
  
 %Estimating current time vector
 n=length(s_i_g_n_a_l);
 t=0:1/sampFreq:(n-1)/sampFreq;
 
 %Estimating new time vector
 tnew = 0:1/newSamplingFreq:(n-1)/sampFreq;
 newTimeSeries = spline(t,s_i_g_n_a_l,tnew);
 
 %exporting new time-series and time vector to MatLab workspace
 assignin('base','t_resampled',tnew);
 assignin('base','sig_resampled',newTimeSeries);
 
end;

% --------------------------------------------------------------------
function varargout = AlignMUAP_Callback(h, eventdata, handles, varargin)


global s_i_g_n_a_l;

%current sampling frequency 
sampFreq = str2double(get(handles.editSampFreq,'String'));

%showing the command window variables in the list box
vars = evalin('base','who');

ss_wnd = evalin('base','ss_wnd');
ss_wpeak = evalin('base','ss_wpeak');

if(isempty(ss_wnd)==1 | isempty(ss_wpeak)==1),
     msgbox('ss_wnd and ss_wpeak are not defined. MUAP peak detection should be performed.','guiDSP','error');
     return;
end;

%Getting aligned MUAPs
[RAW_MUAP,RAW_ALIGNED_MUAP,nMUAPs]= ss_get_alignedMUAPs(ss_wnd,ss_wpeak,s_i_g_n_a_l,sampFreq);

figure(1);
title('Raw MUAPs');

hold on;
for i=1:1:nMUAPs,
    plot(RAW_MUAP{i});
end;
    
figure(2);
title('Aligned MUAPs');
hold on;
for i=1:1:nMUAPs,
    plot(RAW_ALIGNED_MUAP(i,:));
end;

%exporting new time-series and time vector to MatLab workspace
assignin('base','RAW_MUAP',RAW_MUAP);
assignin('base','RAW_ALIGNED_MUAP',RAW_ALIGNED_MUAP);

% --------------------------------------------------------------------
function varargout = IA_MUAPs_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;

%Estimating complex envelope of the input time-series
IA = abs(hilbert(s_i_g_n_a_l));

%current sampling frequency 
sampFreq = str2double(get(handles.editSampFreq,'String'));

%showing the command window variables in the list box
vars = evalin('base','who');

ss_wnd = evalin('base','ss_wnd');
ss_wpeak = evalin('base','ss_wpeak');

if(isempty(ss_wnd)==1 | isempty(ss_wpeak)==1),
     msgbox('ss_wnd and ss_wpeak are not defined. MUAP peak detection should be performed.','guiDSP','error');
     return;
end;

%Getting aligned MUAPs
[IA_MUAP,IA_ALIGNED_MUAP,nMUAPs_IA]= ss_get_alignedMUAPs(ss_wnd,ss_wpeak,IA,sampFreq);

figure(1);
title('Instantaneous amplitude of raw MUAPs');

hold on;
for i=1:1:nMUAPs_IA,
    plot(IA_MUAP{i});
end;
    
figure(2);
title('Instantaneous amplitude of aligned MUAPs');
hold on;
for i=1:1:nMUAPs_IA,
    plot(IA_ALIGNED_MUAP(i,:));
end;

%exporting new time-series and time vector to MatLab workspace
assignin('base','IA_MUAP',IA_MUAP);
assignin('base','IA_ALIGNED_MUAP',IA_ALIGNED_MUAP);

% --------------------------------------------------------------------
function varargout = AR_rawMUAPs_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;

%current sampling frequency 
sampFreq = str2double(get(handles.editSampFreq,'String'));

%showing the command window variables in the list box
vars = evalin('base','who');

ss_wnd = evalin('base','ss_wnd');
ss_wpeak = evalin('base','ss_wpeak');

if(isempty(ss_wnd)==1 | isempty(ss_wpeak)==1),
     msgbox('ss_wnd and ss_wpeak are not defined. MUAP peak detection should be performed.','guiDSP','error');
     return;
end;

dlg_title = 'AR model based on the LMS algorithm';
prompt = {'Order'};
num_lines= 1;
def = {'4'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    order = str2num(msgbox_output{1});
end;

%Estimating AR coefficients
[a,nMUAPs]= ss_get_ARcoefMUAPs(ss_wnd,ss_wpeak,s_i_g_n_a_l,order,0.1,0.001);
%exporting new time-series and time vector to MatLab workspace
assignin('base','AR_RAW_MUAPs',a);
assignin('base','nMUAPs',nMUAPs);



% --------------------------------------------------------------------
function varargout = AR_IA_MUAPs_Callback(h, eventdata, handles, varargin)

global s_i_g_n_a_l;

%Estimating complex envelope of the input time-series
IA = abs(hilbert(s_i_g_n_a_l));

%current sampling frequency 
sampFreq = str2double(get(handles.editSampFreq,'String'));

%showing the command window variables in the list box
vars = evalin('base','who');

ss_wnd = evalin('base','ss_wnd');
ss_wpeak = evalin('base','ss_wpeak');

if(isempty(ss_wnd)==1 | isempty(ss_wpeak)==1),
     msgbox('ss_wnd and ss_wpeak are not defined. MUAP peak detection should be performed.','guiDSP','error');
     return;
end;

[a,nMUAPs]= ss_get_ARcoefMUAPs(ss_wnd,ss_wpeak,IA,4,0.1,0.001);
%exporting new time-series and time vector to MatLab workspace
assignin('base','AR_IA_MUAPs',a);
assignin('base','nMUAPs',nMUAPs);


% --------------------------------------------------------------------
function Butterworth_Callback(hObject, eventdata, handles)
% hObject    handle to Butterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function HPButterworth_Callback(hObject, eventdata, handles)
% hObject    handle to HPButterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s_i_g_n_a_l;

prompt = {'Cutoff frequency (Hz)', 'order'};
dlg_title = 'High-pass Butterworth filter';
num_lines= 1;
def = {'20','4'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    cutoffFreq = str2num(msgbox_output{1});
    order = str2num(msgbox_output{2});
    %current sampling frequency 
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    Wn = cutoffFreq/(sampFreq/2);
    [b,a] = butter(order,Wn,'high'); %filter design
    figure(1); freqz(b,a,sampFreq/2,sampFreq); %plotting filter characteristics
    y = filtfilt(b,a,s_i_g_n_a_l); %zero-phase digital filtering (zero-phase distortion)
    %exporting filtered signal
    assignin('base','HighPassButFiltSig',y);
end;
% --------------------------------------------------------------------
function PassBandButterworth_Callback(hObject, eventdata, handles)
% hObject    handle to PassBandButterworth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global s_i_g_n_a_l;

prompt = {'Lower cutoff frequency (Hz)', 'Upper cutoff frequency (Hz)','order'};
dlg_title = 'Pass-band Butterworth filter';
num_lines= 1;
def = {'20','1500','4'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    LcutoffFreq = str2num(msgbox_output{1});
    UcutoffFreq = str2num(msgbox_output{2});
    order = str2num(msgbox_output{3});
    %current sampling frequency 
    sampFreq = str2double(get(handles.editSampFreq,'String'));
    Wn = [LcutoffFreq UcutoffFreq]/(sampFreq/2);
    [b,a] = butter(order,Wn); %filter design
    figure(1); freqz(b,a,sampFreq/2,sampFreq); %plotting filter characteristics
    y = filtfilt(b,a,s_i_g_n_a_l); %zero-phase digital filtering (zero-phase distortion)
    %exporting filtered signal
    assignin('base','PassBandButFiltSig',y);
end;

% --------------------------------------------------------------------
function varargout = ss_DigTrigger_Callback(h, eventdata, handles, varargin)

DigitalTrigger;%calling window


% --------------------------------------------------------------------
function CleanBackActivity_Callback(hObject, eventdata, handles)
% hObject    handle to CleanBackActivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global s_i_g_n_a_l;

ss_wnd = evalin('base','ss_wnd');
ss_wpeak = evalin('base','ss_wpeak');

if(isempty(ss_wnd)==1 | isempty(ss_wpeak)==1),
     msgbox('ss_wnd and ss_wpeak are not defined. MUAP peak detection should be performed.','guiDSP','error');
     return;
end;

prompt = {'Threshold - peak-to-peak (V)'};
dlg_title = 'Cleaning background activity';
num_lines= 1;
def = {'0.2'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
threshold = str2num(msgbox_output{1});
end;

%Eliminating background activity
[wpeak_updated,wnd_update]= ss_backgroundActRemoval(ss_wpeak,ss_wnd, s_i_g_n_a_l,threshold);

%exporting new time-series and time vector to MatLab workspace
assignin('base','ss_wnd',wnd_update);
assignin('base','ss_wpeak',wpeak_updated);


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_7_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function EstimatePCs_Callback(hObject, eventdata, handles)
% hObject    handle to EstimatePCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PlotPCs_Callback(hObject, eventdata, handles)
% hObject    handle to PlotPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function VarPCs_Callback(hObject, eventdata, handles)
% hObject    handle to VarPCs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function PlotScores_Callback(hObject, eventdata, handles)
% hObject    handle to PlotScores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------
function PCAanalysis_Callback(hObject, eventdata, handles)

PCAanalysis; %calling dialog


% --------------------------------------------------------------------
function MUAPvisualization_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
if(isempty(MU)==1),
     msgbox('MU is not defined.','guiDSP','error');
     return;
end;

VisualizeMUAPtemplates(MU);

% --------------------------------------------------------------------
function MUAP_sim_Callback(hObject, eventdata, handles)

figure(1987); %creating figure
MU = zeros(1,1); %initializing variable
[MU,ElecPos] = SimulateMU;
%exporting motor unit structure to workspace
assignin('base','MU',MU);
assignin('base','ElecPos',ElecPos);

% --------------------------------------------------------------------
function MUAP_Pat_Gen_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
if(isempty(MU)==1),
     msgbox('MU is not defined.','guiDSP','error');
     return;
end;

[MUAPpatterns] = GenerateMUAPpatterns(MU);
%exporting motor unit structure to workspace
assignin('base','MUAPpatterns',MUAPpatterns);


% --------------------------------------------------------------------
function MUAPT_sim_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
if(isempty(MU)==1),
     msgbox('MU is not defined.','guiDSP','error');
     return;
end;

[MU,t,Fs] = GenerateMUAPT(MU);

%exporting variables to workspace
assignin('base','MU',MU);
assignin('base','t',t);
assignin('base','Fs',Fs);

% --------------------------------------------------------------------
function MUAPT_visualization_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
t = evalin('base','t');
Fs = evalin('base','Fs');

if(isempty(MU)==1 | isempty(t)==1 | isempty(Fs)==1  ),
     msgbox('MU/t/Fs is not defined.','guiDSP','error');
     return;
end;

VisualizeMUAPT(MU,t,Fs);

% --------------------------------------------------------------------
function Firing_time_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
t = evalin('base','t');
Fs = evalin('base','Fs');

if(isempty(MU)==1 | isempty(t)==1 | isempty(Fs)==1  ),
     msgbox('MU/t/Fs is not defined.','guiDSP','error');
     return;
end;

VisualizeMUAPTfiringTime(MU,t,Fs);
% --------------------------------------------------------------------
function IPI_histogram_Callback(hObject, eventdata, handles)

prompt = {'Number of bins'};
dlg_title = 'Inter-spike interval histogram';
num_lines= 1;
def = {'10'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    Nbins = str2num(msgbox_output{1});
end;


%Getting motor unit structure
MU = evalin('base','MU');

if(isempty(MU)==1  ),
     msgbox('MU is not defined.','guiDSP','error');
     return;
end;

MUAPT_IPIhist(MU,Nbins);

% --------------------------------------------------------------------
function MuscleCrossSec_Callback(hObject, eventdata, handles)

MU = evalin('base','MU');
ElecPos = evalin('base','ElecPos');

if(isempty(MU)==1 | isempty(ElecPos)==1),
     msgbox('MU/ElecPos is not defined.','guiDSP','error');
     return;
end;
figure(111);
VisualizeMuscle(MU,ElecPos);


% --------------------------------------------------------------------


% --------------------------------------------------------------------
function MUAPOverlapBld_Callback(hObject, eventdata, handles)
% hObject    handle to MUAPOverlapBld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%calling GUI
MUAPOverlapBuilder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GenMUAPT_Callback(hObject, eventdata, handles)

%Creating input dialog box
dlg_title = 'MUAPT simulation';
prompt = {'t (ms)','Fs (Hz)','Number of MUs'};
num_lines= 1;
def = {'1000','10040','5'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    t = str2num(msgbox_output{1}); %time (ms)
    Fs = str2num(msgbox_output{2}); %sampling frequency (Hz)
    Nmus = str2num(msgbox_output{3});  %number of motor units
end;

%generating motor unit action potential trains
[MU,t,Fs] = exp_GenerateMUAPT(t,Fs,Nmus);
%exporting variables to workspace
assignin('base','MU',MU);
assignin('base','t',t);
assignin('base','Fs',Fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GenSEMG_Callback(hObject, eventdata, handles)

%Creating input dialog box
dlg_title = 'sEMG simulation';
prompt = {'t (ms)','Fs (Hz)','Number of MUs','SNR (dB)','Gain of Amplifier'};
num_lines= 1;
def = {'1000','10040','5','20','10000'};
msgbox_output = inputdlg(prompt,dlg_title,num_lines,def);

if(isempty(msgbox_output)==1),
    return;
else
    t = str2num(msgbox_output{1}); %time (ms)
    Fs = str2num(msgbox_output{2}); %sampling frequency (Hz)
    Nmus = str2num(msgbox_output{3});  %number of motor units
    SNR = str2num(msgbox_output{4});  %signal-to-noise ration (dB)
    AmplifierGain = str2num(msgbox_output{5}); %Gain of amplifier
end;

%generating motor unit action potential trains
[MU,t,Fs] = exp_GenerateMUAPT(t,Fs,Nmus);

%simulating surface EMG signal
[EMG] = exp_GenerateEMG(MU,SNR,AmplifierGain);

%exporting variables to workspace
assignin('base','MU',MU);
assignin('base','t',t);
assignin('base','Fs',Fs);
assignin('base','EMG',EMG);
assignin('base','SNR',SNR);
assignin('base','AmplifierGain',AmplifierGain);
