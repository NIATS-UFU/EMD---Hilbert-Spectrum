% Function name....: extrema
% Date.............: December 02, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Extrema estimates local maxima (peaks) and minima (valleys) of time-series.
% Parameters.......: 
%                    Y -> input time-series
% Return...........:
%                    LocalMaxV -> Local maxima
%                    LocalMaxI -> Local maximum indexes
%                    LocalMinV -> Local minima
%                    LocalMinI -> Local minimum indexes
% mcc -x extrema

function [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(Y)

Len= max(size(Y)); % time-series length 
    
EXTREMA = zeros(4,Len);
    
LocalMaxV=[];        
LocalMaxI=[];        
LocalMinV=[];
LocalMinI=[];

counter_1 = 0;
counter_2 = 0;
    
for x=2:Len-1
    if Y(x)>Y(x-1) & Y(x)>=Y(x+1)
        counter_1 = counter_1 + 1;
%       LocalMaxI=[LocalMaxI; x];     % index
%       LocalMaxV=[LocalMaxV; Y(x)];  %value
        EXTREMA(1,counter_1) = x;
        EXTREMA(2,counter_1) = Y(x);
    elseif  Y(x)<Y(x-1) & Y(x)<=Y(x+1)
        counter_2 = counter_2 + 1;
%       LocalMinI=[LocalMinI; x];      % index
%       LocalMinV=[LocalMinV; Y(x)];  % value
        EXTREMA(3,counter_2) = x;
        EXTREMA(4,counter_2) = Y(x);
       end;
end; %function 

LocalMaxI = EXTREMA(1,1:counter_1)';
LocalMaxV = EXTREMA(2,1:counter_1)';

LocalMinI = EXTREMA(3,1:counter_2)';
LocalMinV = EXTREMA(4,1:counter_2)';