% Function name....: meanEnv
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    meanEnv estimates the average envelope of a time-series. 
% Parameters.......: 
%                    Y ......... -> input time-series
% Return...........:
%                    AveEnvelope -> average envelope
%                    yu ........ -> upper envelope
%                    yl ........ -> lower envelope
%                    n_extrema.. -> the sum of the number of peaks and valleys
%
% Remarks..........:
%                    meanEnv uses the extrema function to estimate the local maxima and minima
%                    of the signal Y and after it uses a cubic spline interpolation to
%                    connect those points in order to obtain the lower and upper envelope.
%                    With the aim of minimizing the end effect in the spline fitting a 
%                    solution is provided.

function [env]= EnvMAF(y,pulseSize)

    sizeY = length(y);
    
    %Generating a rectangular pulse with area equals 1
    t=1:1:sizeY;
    rectPulse = rectpuls(t,2*pulseSize)/pulseSize;
    
    %Calculating the DFT of the rectangular pulse
    rectDFT = fft(rectPulse);
    
    %Rectifying the input signal y
    yr = abs(y);
    
    %Calculating the DFT of the rectified signal
    yrDFT = fft(yr);
    
    %Filtering the input signal
    yfilt = rectDFT.*yrDFT';
    
    %Estimating the IFFT of the filtered signal
    env = abs(ifft(yfilt));
    
    
 