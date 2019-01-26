% Function name....: meanEnv
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    meanEnv estimates the average envelope of a time-series. 
% Parameters.......: 
%                    Y ............. -> input time-series
%                    interpMethod .. -> interpolation method (0=spline, ~=0 pchip)
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
% For compilation > mcc -x meanEnv
function [AveEnvelope,yu,yl,n_extrema]= meanEnv(Y,interpMethod)

    %initializing vectors with null
    AveEnvelope=[];
    yu=[];
    yl=[];
    n_extrema = 0;
    
    %estimating upper and lower envelope of the time-series Y
    [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(Y);%estimating signal upper and lower extrema
    if (isempty(LocalMaxV)|isempty(LocalMinV)),
        return;
    end%if

    %Calculating the number of extrema
    n_extrema=length(LocalMaxV)+length(LocalMinV);
    
    %Performing end corrections for the upper envelope
    s=length(Y);
    x1=LocalMaxV(1);	
    x2=LocalMaxV(length(LocalMaxV));
   
    % The extreme points of the averaged envelope are points from the original sequence
    xx=[x1 LocalMaxV' x2]; 
    t=1:s; %buiding time vector
    tt=[t(1)-1/300 t t(s)+1/300];%rearranging time vector in order to have additional samples
    ind=[t(1)-1/300 LocalMaxI' t(s)+1/300];
    
    if(interpMethod~=0),
        yu=pchip(ind,xx,tt);
    else
        yu=spline(ind,xx,tt);
    end;
    yu = yu(2:s+1);
    
    % Performing end corrections for the lower envelope
    x1=LocalMinV(1);	
    x2=LocalMinV(length(LocalMinV));
    xx=[x1 LocalMinV' x2]; 
    ind=[t(1)-1/300 LocalMinI' t(s)+1/300];
    
    if(interpMethod~=0),
         yl=pchip(ind,xx,tt);
    else
         yl=spline(ind,xx,tt);
    end;
    yl = yl(2:s+1);
    
    %Estimating the average envelope
    AveEnvelope = (yu+yl)/2;
 