% Function name....: FresponseWLPD (weighted low pass differential filter)
% Date.............: October 09, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Frequency response of a low-pass differential filter. 
% Parameters.......: 
%                    f ............. -> frequency (vector)
%                    fs............. -> sampling frequency (Hz)
%                    N.............. -> window size 
% Return...........:
%                    H -> Gain as function of f

function [H] = FresponseWLPD(f,N,fsr)

nsamples=length(f);

for i=1:fsr/2,
    sum=0;
    for nn=1:N,
        sum= sum + ((exp(2j*pi*f(i)/fsr*nn)) - (exp(-2j*pi*f(i)/fsr*nn))) ;
    end%for nn
    H(i) = abs(sum)/N;
end;%for
