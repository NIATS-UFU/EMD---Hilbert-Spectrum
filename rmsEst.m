% Function name....: rmsEst
% Date.............: September 11, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    rmsEst estimates the root mean square (RMS) of an input vector and return the central
%                    position of the window used in the estimation.
%                    
% Parameters.......: 
%                    vector .........-> input vector
%                    windowSize .....-> size of the window for RMS estimation
% Return...........:
%                    rms....-> root-mean-square value of a sequence of windows of size windowSize
%                    index..-> a sequence of indexes indicating the position of  estimated RMS values

function [rms,index] = rmsEst(vector,windowSize)

VectorSize = length(vector);
N = fix(VectorSize/windowSize);

h = waitbar(0,'Please wait...');
for i=1:N,
   waitbar(i/N);
   rms(i) =norm(vector(windowSize*(i-1)+1:i*windowSize))/sqrt(windowSize);
   index(i) = round(windowSize*(i-1) + windowSize/2);
end;
close(h);%closing waiting bar