function [PassBandButFiltSig,i_m_fs]= EstIMFsBatch(Y)
% Filtering
Wn = [20 1000]/(10040/2);
[b,a] = butter(4,Wn); %filter design
PassBandButFiltSig = filtfilt(b,a,Y); %zero-phase digital filtering (zero-phase distortion)

%Estimating IMFs
[i_m_fs,residue]= sig_to_imf(PassBandButFiltSig,1e-5,1);