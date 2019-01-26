% Function name....: meanEnv
% Date.............: December 12, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    meanEnv estimates the average envelope of a time-series. 
% Parameters.......: 
%                    Y .....-> input time-series
%                    r .....-> resample factor (it must be an integer greater than one)
%                    method -> It should be one of the following interpolation methods:
%                               * 'nearest': Nearest neighbor interpolation
%                               * 'linear'.: Linear interpolation 
%                               * 'spline'.: Cubic spline interpolation
%                               * 'pchip'..: Piecewise cubic Hermite interpolation
% Return...........:
%                    meanAmpVecEnvelope  -> average envelope
%                    UpperEnvelope       -> upper envelope
%                    LowerEnvelope       -> lower envelope
%                    UpperEnv_index      -> upper envelope indexes (e.g. plot(UpperEnv_index,UpperEnvelope))
%                    LowerEnv_index      -> lower envelope indexes (e.g. plot(LowerEnv_index,LowerEnvelope))
%
% Remarks..........:
%                    meanEnv uses the envelope function to estimate the upper and lower envelope of the
%                    input signal. After calculating the average envelope it is resampled via the cubic
%                    spline interpolotion in order to have the same number of samples of the Y signal.

function [meanAmpVecEnvelope,...
         UpperEnvelope,UpperEnv_index,LowerEnvelope,LowerEnv_index]= meanEnv(Y,r,method)

   %initializing vectors with null
    meanAmpVecEnvelope=[];
    meanAmpVecEnvelope_index=[];
    UpperEnvelope=[];
    UpperEnv_index=[];
    LowerEnvelope=[];
    LowerEnv_index=[];
    
    %estimating upper and lower envelope of the time-series Y
    [UpperEnvelope,UpperEnv_index,LowerEnvelope,LowerEnv_index]= envelope(Y,r,method);
    
    %determining the size of each envelope vector
    vecSizeUpper = max(size(UpperEnvelope));
    vecSizeLower = max(size(LowerEnvelope));
    
    %defining reference and auxiliar vectors 
    %the reference vector is the one that has less number of samples
    if (vecSizeUpper > vecSizeLower) 
        referenceVec = LowerEnv_index;
        auxVec = UpperEnv_index;
        referenceAmpVec = LowerEnvelope;
        auxAmpVec = UpperEnvelope;
    else
        referenceVec = UpperEnv_index;
        auxVec = LowerEnv_index;
        referenceAmpVec = UpperEnvelope ;
        auxAmpVec = LowerEnvelope;        
    end %if
    
    % defining the size of the reference and auxiliary vectors
    sizeReferenceVec = max(size(referenceVec));
    sizeAuxVec = max(size(auxVec));
    meanAmpVecEnvelope_index=referenceVec;
    
    %initializing auxiliary variables used in the averaging process 
    %It is important to notice that two sequences with different numbers of samples
    %are being added. This is an aproximation that works very well.
    p1 = auxVec(1);
    p2 = auxVec(2);
    k = 3;
    
    for I = 1:sizeReferenceVec,
        if (referenceVec(I)<p1),
                meanAmpVecEnvelope(I)= referenceAmpVec(I)/2;
        elseif (referenceVec(I)>p2) & (k >=sizeAuxVec),
                meanAmpVecEnvelope(I)= referenceAmpVec(I)/2;
        else
                while not ((p1 <= referenceVec(I)) & (p2>=referenceVec(I))),
                      p1 = p2;
                      p2 = auxVec(k);
                      k=k+1;
                      if(k >=sizeAuxVec) 
                          break;
                      end%if
                end %while
                if(abs(p1-referenceVec(I)) <= abs(p2-referenceVec(I))),
                    meanAmpVecEnvelope(I)= (auxAmpVec(k-2) + referenceAmpVec(I))/2;
                else
                    meanAmpVecEnvelope(I)= (auxAmpVec(k-1) + referenceAmpVec(I))/2;
                end %if
        end%switch
   end%for
      
   %The final step is to perform end corrections for splines over the averaged envelope
   s=length(Y);
   x1=Y(1);	
   x2=Y(s);
   % The extreme points of the averaged envelope are points from the original sequence
   xx=[x1 meanAmpVecEnvelope x2]; 
   t=1:s; %buiding time vector
   tt=[t(1)-1/300 t t(s)+1/300];%rearranging time vector in order to have to additional samples
   ind=[t(1)-1/300 meanAmpVecEnvelope_index t(s)+1/300];
   ys=spline(ind,xx,tt);
   meanAmpVecEnvelope = ys(2:s+1);

  