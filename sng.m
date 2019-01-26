% Function name....: sng
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    sng is defined as:
%
%                    sng(x)= 1, if x>0
%                    sng(x)= 0, otherwise       
%
% Parameters.......: 
%                    x .....-> a scalar input 
% Return...........:
%                    s  -> outuput as defined above
% for compilation: mcc -x sng

function [s]= sng(x)

  if x>0,
      s=1;
  else
      s=0;
  end;

   