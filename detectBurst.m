% Function name....: threshold
% Date.............: January 02, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Envelope etimates upper and lower envelopes of time-series. 
% Parameters.......: 
%                    Y .....-> input time-series
%                    r .....-> resample factor (it must be an integer greater than one)
%                    method -> It should be one of the following interpolation methods:
%                               * 'nearest': Nearest neighbor interpolation
%                               * 'linear'.: Linear interpolation 
%                               * 'spline'.: Cubic spline interpolation
%                               * 'pchip'..: Piecewise cubic Hermite interpolation
% Return...........:
%                    UpperEnvelope  -> upper envelope
%                    LowerEnvelope  -> lower envelope
%                    UpperEnv_index -> upper envelope indexes (e.g. plot(UpperEnv_index,UpperEnvelope))
%                    LowerEnv_index -> lower envelope indexes (e.g. plot(LowerEnv_index,LowerEnvelope))
%

function [onset,tr]= detectBurst(y,MaxRange,Nbits,fs)

onset = [];

[meanG,index_meanG,sig,sigIndex,tr]= threshold(y,MaxRange,Nbits);

ySize = length(y);
yrect = abs(y);

yMax = max(y);
yMin = min(y);

for i=1:ySize,

    if (yrect(i)>=tr)
        onset(i) = yMax;
    else
        onset(i) = yMin;
    end%if
    
end%for