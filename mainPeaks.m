% Function name....: mainPeaks
% Date.............: July 28, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This function estimates the principal peaks of an input time-series
% Parameters.......: 
%                    s_i_g_n_a_l .....-> input time-series
%                    maxADRange ......-> maximum AD range
%                    Nbits............-> number of bits (resolution)
%                    sampFreq.........-> sample frequency (Hz)
% Return...........:
%                    Index .... -> localization in time when principal peaks occur
%                    Amplitude  -> amplitude of principal peaks at a especific time

% Remarks: detectBurst and extrema are functions already implemented and do not belong to the MatLab function library.

function [Index,Amplitude]= mainPeaks(s_i_g_n_a_l,maxADRange,Nbits,sampFreq)

    Index = [];
    Amplitude = [];
    RelevantLocalMaxIndex = [];
    
    %Estimating signal offset
    [onset,tr]= detectBurst(s_i_g_n_a_l,maxADRange,Nbits,sampFreq);

    %calculating signal extrema (valleys and peaks).
    [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(s_i_g_n_a_l);
    
    %Estimating principal peaks
    RelevantLocalMaxIndex = find(LocalMaxV>=tr); 
    Index = LocalMaxI(RelevantLocalMaxIndex);
    Amplitude = LocalMaxV(RelevantLocalMaxIndex);