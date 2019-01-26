% Function name....: envelope
% Date.............: December 03, 2002
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

function [UpperEnvelope,UpperEnv_index,LowerEnvelope,LowerEnv_index]= envelope(Y,r,method)

    
    UpperEnvelope=[];  % initiating upper envelope vector = null
    LowerEnvelope=[];  % initiating lower envelope vector = null
    UpperEnv_index=[];
    LowerEnv_index=[];
    xi_upper=[];  
    xi_lower=[];
    
    [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(Y);%estimating signal upper and lower extrema
    
    %%%%%%%%%%%%%%%%%%
    LenMax= max(size(LocalMaxV)); % time-series length 
    xi_upper=[1:1/r:LenMax]; %generating resampled index vector
    
    LenMin= max(size(LocalMinV)); % time-series length 
    xi_lower=[1:1/r:LenMin]; %generating resampled index vector
        
    UpperEnvelope = interp1(LocalMaxV,xi_upper,method); %estimating upper envelope
    LowerEnvelope = interp1(LocalMinV,xi_lower,method); %estimating lower envelope

    %%%%%%%%%%%%%%%%%%
    UpperEnv_index = interp1(LocalMaxI,xi_upper,method); %estimating upper envelope indexes 
    LowerEnv_index = interp1(LocalMinI,xi_lower,method); % estimating lower envelope indexes