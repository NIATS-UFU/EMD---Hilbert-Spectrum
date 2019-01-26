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

function [a]= meanEnvelope(Y,r,method)

   % meanAmpVecEnvelope=[];
    
    UpperEnvelope=[];
    UpperEnv_index=[];
    LowerEnvelope=[];
    LowerEnv_index=[];
    
    
    %estimating upper and lower envelope of the time-series Y
 %   [UpperEnvelope,UpperEnv_index,LowerEnvelope,LowerEnv_index]= envelope(Y,r,'cubic');
    
    %determining the size of each envelope vector
 %   vecSizeUpper = size(max(UpperEnvelope));
 %   vecSizeLower = size(max(LowerEnvelope));
    
    %defining reference and auxiliar vectors 
 %   if (vecSizeUpper > vecSizeLower) 
 %       referenceVec = LowerEnv_index;
 %       auxVec = UpperEnv_index;
 %       referenceAmpVec = LowerEnvelope;
 %       auxAmpVec = UpperEnvelope;
 %   else
 %       referenceVec = UpperEnv_index;
 %       auxVec = LowerEnv_index;
 %       referenceAmpVec = UpperEnvelope ;
 %       auxAmpVec = LowerEnvelope;        
 %   end %if
    
    % defining the size of the reference vector
 %   sizeReferenceVec = size(max(referenceVec));
    
 %   p1 = auxVec(1);
 %   p2 = auxVec(2);
 %   k = 3;
    
 %   for I = 1:sizeReferenceVec,
        
 %       while NOT ((p1 <= referenceVec(I)) & (p2>=referenceVec(I))),
 %           p1 = p2;
 %           p2 = auxVec(k);
 %           k=k+1;
 %       end %while
        
 %       if(abs(p1-referenceVec(I)) <= abs(p2-referenceVec(I))),
 %           meanAmpVecEnvelope(I)= (auxAmpVec(k-2) + referenceAmpVec(I))/2;
 %       else
 %           meanAmpVecEnvelope(I)= (auxAmpVec(k-1) + referenceAmpVec(I))/2;
 %       end %if
 %   end%for




  