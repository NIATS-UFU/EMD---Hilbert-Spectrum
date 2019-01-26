% Function name....: histogram
% Date.............: March 30, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This function estimates the histogram of a time
%                    series. This histogram has Nbits different levels,
%                    which are based on the number of bits of the A/D
%                    converter.
% Parameters.......: 
%                    Y .....-> input time-series
%                Nbits .....-> number of bits of the A/D converter
% Return...........:
%                    xii  -> samples (amplitude)
%                    h  -> number of occurrences (frequency)

function [xii,h]= histogram(y,Nbits)

%Estimating signal histogram
[N,X] = hist(y,2^Nbits);
ii = find(N);
xii = X(ii);
h = N(ii);

%plot(X(ii),N(ii))
% %Estimating the histogram
%  x=MinRange:(MaxRange-MinRange)/2^Nbits:MaxRange;
%  [n,xout] = hist(y,x);
%  
%  %Eliminating zeros and interpolating the histogram in order to increase the precision
%  i=find(n);
%  [LocalMaxV,LocalMaxI, LocalMinV,LocalMinI]= extrema(n(i));
% 
%  ii=i(LocalMaxI);
%  xii= xout(ii(1)):MaxRange/2^Nbits:xout(ii(length(ii)));
%  h = interp1(xout(ii),LocalMaxV,xii,'pchip');
 
