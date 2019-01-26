% Function name....: IPIplot
% Date.............: July 01, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    IPIplot plots the instant when firings of simulated MUs occur
% Parameters.......: 
%                    t ........-> vector containing the time when firings occur
%                    ymin .....-> minimum amplitude 
%                    ymax .....-> max amplitude
%                    color .....-> one of the possible colors available 
% Example:  IPIplot(IPI(1,:),min(emg),max(emg),'r');

function []= IPIplot(t,ymin,ymax,color)
 
    S=max(size(t));
    for i=1:S,
        if(t(i)>0)
           hold on; plot([t(i) t(i)],[ymin ymax],color);
       end;
    end;
    