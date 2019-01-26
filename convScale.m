% Function name....: convScale
% Date.............: December 28, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    convScale is useful for converting value from scale x into y .
% Parameters.......: 
%                    xmin .....-> minimum value of x scale
%                    xmax .....-> maximum value of x scale
%                    x    .....-> value to be converted to y scale
%                    ymin .....-> minimum value of y scale
%                    ymax .....-> maximum value of y scale
% Return...........:
%                    y    .....-> converted value

function [y]= convScale(xmin,xmax,x,ymin,ymax)

     y=(x-xmin)*(ymax-ymin)./(xmax-xmin)+ymin;
     
     
     