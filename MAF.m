% Function name....: MAF
% Date.............: September 11, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    MAF is the implementation of a moving average filter
%                    
% Parameters.......: 
%                    signal .........-> input signal
%                    M......... .....-> filter order
% Return...........:
%                    FilteredSignal....-> filtered signal

function [FilteredSignal]= MAF (signal,M)

    FilteredSignal = [];
 
    nrAmos = length(signal);
    
    h = waitbar(0,'Please wait...');
    for n=1:nrAmos,
        waitbar(n   /nrAmos);
        
        sum = 0;
        for k=1:M,
            if (n-k)>0,
                sum = sum + signal(n-k);
            else
                sum = signal(n)*M;   
                break;
            end%if
        end %for k
    
        FilteredSignal(n) = sum / M;
    
    end%for n
    close(h);%closing waiting bar