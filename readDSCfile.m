% Function name....: readDSCfile
% Date.............: February 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This functions estimates the power spectrum of a time-series based on its AR model 
% Parameters.......: 
%                    a...........-> coefficients of the AR model
%                    bins........-> number of samples for frequency estimation
%                    Fs..........-> sampling frequency
%                    df .........-> frequency resolution (Hz). When using AR models for estimating the power spectrum of signals
%                                   we do not have the dependence on the size of the signal window in analysis, i.e. this
%                                   variable is not definied as df = Fs/N like in Fourier analysis.
% Return...........:
%                    freq........-> frequency in Hz
%                    s...........-> normalized power spectrum
%
% Reference: Akay, Metin, "Biomedical Signal Processing", Academic Press, pp. 204 - 206, 1994
% Note: for compilation use 'mcc -x readDSCfile'

function [forceLevel]=readDSCfile(filename,forceSig,Fs)

Nsamples = length(forceSig);

FID = fopen(filename,'rt'); %opening the .dsc file

f = zeros(1,10);


TLINE = 1;

%Reading force level
k=0;
if FID>0,
    for i=1:1:10,
        TLINE = fgets(FID);
        if(i>5),
            k = k+1;
            [token,rem] = strtok(TLINE,':');
            [token,rem] = strtok(rem,':');
            [token,rem] = strtok(token,'to');
            f(k) = str2num(token);
            [token,rem] = strtok(rem,'to');
            k = k+1;
            f(k) = str2num(token);
        end;
    end;
    fclose(FID); %closing file
else
    disp('cannot open file');
    forceLevel = zeros(1,1);
    return;
end;
disp(f);


w = round(0.01 * Fs); %number of samples
forceLevel = zeros(1,Nsamples-w);


for n=1:1:Nsamples-w,
    
   % disp(n);
    meanF = mean(forceSig(n:n+w));
    
    if (meanF>=f(1)) & (meanF<f(2)),
        forceLevel(n) = 1; %first level
    end;
    
    if (meanF>=f(3)) & (meanF<f(4)),
        forceLevel(n) = 2; %second level
    end;   
    
    if (meanF>=f(5)) & (meanF<f(6)),
        forceLevel(n) = 3; %third level
    end;

    if (meanF>=f(7)) & (meanF<f(8)),
        forceLevel(n) = 4; %fourth level
    end;    
    
    if (meanF>=f(9)) & (meanF<=f(10)),
        forceLevel(n) = 5; %fith level
    end;    
    
    if (meanF<f(1)),
        forceLevel(n) = 0; %bellow the first level
    end;
    
    if (meanF>f(10)),
        forceLevel(n) = 6; %above the fith level
    end;        
end;%for