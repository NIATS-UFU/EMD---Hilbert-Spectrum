% Function name....: WLPD (weighted low pass differential filter)
% Date.............: April 24, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This filter is designed to accentuate the MUAP spike in SEMG signals and
%                    supress background activity. It is based on the work 'Digital filter design 
%                    for peak detection of surface EMG', Zhengquan Xu & Shaojun Xiao,
%                    Journal of Electromyography and Kinesiology, 10 (2000) 275-281
% Parameters.......: 
%                    x ............. -> input time-series
%                    fs............. -> sampling frequency
%                    N.............. -> window size 
% Return...........:
%                    y -> filtered signal

function [y] = WLPD(x,fs,N)

y = [];
nsamples = length(x); 
Ns = round(N*fs/10000); %window size (details can be found in the work of Xu/Xiao)
%estimating sinusoidal window
k=1:1:Ns;
ws = sin(pi*k./Ns);

h = waitbar(0,'Please wait...');
%Filtering
for i=1:nsamples,
    waitbar(i/nsamples);
    sum=0;
    for j=1:Ns,
        if(i+j)>nsamples,
            aux1 = x(nsamples);
        else
            aux1 = x(i+j);
        end;%if
        if(i-j)<1,
            aux2 = x(1);
       %     aux2 = 0;
        else
%             aux2 = 0;
           aux2 = x(i-j);
            
        end; %if
       sum = sum + ws(j)*(aux1-aux2);
  %     sum = sum + (aux1-aux2);

    end;%for
    y(i)=sum; %filtered signal
end;%for
close(h);%closing waiting bar


