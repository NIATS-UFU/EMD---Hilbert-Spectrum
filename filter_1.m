
% signal: input signal
% fsr: sampling frequency
% fc: cut-off frequency
% n: order

function [filterCoef]= filter_1 (signal,fsr,fc,n)

filterCoef = [];

nrAmos = lenght(signal);

%Filter coeficients 
for f=0:nrAmos-1 ,
    % f and fsr (Hz)
	w = (2*pi*f)/fsr; %(rad/s)
	filterCoef[f] = sqrt(1/(1+(((1/(tan(pi*fc/fsr)))*tan(w/2))^(2*n))));
end