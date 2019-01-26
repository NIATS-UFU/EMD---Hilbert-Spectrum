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
% for compilation: mcc -x zero_crossing

function [z]= zero_crossing(y)

  z=0;
  s=length(y);
  for I=1:s-1,
      
      %This piece of code replaced the function sng
      x = -y(I)*y(I+1);
      if(x>0),
        z=z+1;
      end;
%     z=z+sng(-y(I)*y(I+1));
      if(y(I+1)==0)&(y(I)~=0),
          z=z+1; %consider that case of sawtooth function (ie, sawtooth(2*pi*50*t))
      end;
  end%for
  
