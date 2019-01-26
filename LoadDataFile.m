% Function name....: ss_get_ARcoefMUAPs
% Date.............: March 08, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Extracted AR coefficients from MUAPs 
%
% Parameters.......: 
%                    ss_wnd ..... -> definition of the boundaries that contain a candidate MUAP
%                    ss_wpeak.... -> peaks and their position in time
%                    signal...... -> input time-series
%                    order....... -> order of the autoregressive model
%                    stepSize.... -> constant of convergence of the LMS algorithm
%                    tolerance... -> minimum mean squared error
% Return...........:
%                    a............-> autoregressive coefficients
%                    nMUAPs............-> number of detected MUAPs
% %Note: for compilation use 'mcc -x LoadDataFile'

function [DATA]= LoadDataFile(fileName)

%opening file
fid = fopen(fileName,'rt');

if (fid<0),
    fprintf('Cannot open file');
    return;
end;

i = 0;
while 1,
    disp(i);
    S = fgetl(fid);
    if(S == -1),
        return;
    end;

    N = str2num(S);
    
    if (isempty(N)~=1),
        i = i+1;
        DATA(i,:) = N;
    end;
    
end; %while