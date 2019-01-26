% Function name....: instantAtrib
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    instantAtrib is useful in calculating instantaneous atributes of a time series.
% Parameters.......: 
%                    signal .....-> input time-series 
%                    fs     .....-> sampling frequency in Hz
% Return...........:
%                    amp    .....-> instantaneous amplitude
%                    phase  .....-> instantaneous phase
%                    freq   .....-> instantaneous frequency

function [amp,phase,freq]= instantAtrib(signal,fs)

     dt=1/fs;%time resolution
     
     y=hilbert(signal);%generating the analytical signal
     
     amp = abs(y); %instantaneous amplitude
     
     %it is important to unwrape the phase in order to eliminate discontinuaties 
     phase = unwrap(angle(y)); %instantaneous phase (in radians)
     
     %padding the phase sequence to obtain the instantaneous frequency
     aux = [phase phase(length(phase))];
     freq=(abs(diff(aux))/(2*pi*dt)); %frequency in Hertz
     freq(length(freq))=freq(length(freq)-1);