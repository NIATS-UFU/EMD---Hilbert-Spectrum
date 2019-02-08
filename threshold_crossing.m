% Function name....: zero_crossing
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    zero_crossing estimates the number of zero-crossings of an input sequence y
%
% Parameters.......: 
%                    y ..... -> input sequence
% Return...........:
%                    z ..... -> number of zero-crossings
% Remarks..........: This function uses sng to define a zero crossing

function [z,index]= threshold_crossing(y,threshold)

  z=0;
  k = 0;
  s=length(y);
  for I=1:s-1,
      if(y(I+1)>=threshold)&&(y(I)<threshold),%---> - to + (direction)
          z=z+1; 
          k = k+1;
          index(k) = I;
      end;
      if(y(I+1)<threshold)&&(y(I)>=threshold),%<--- + to - (direction)
          z=z+1; 
          k = k+1;
          index(k) = I;
      end;      
  end%for
  

   