% Function name....: MHSpec
% Date.............: December 28, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    plotHS generates the Hilbert Spectrum of a signal.
% Parameters.......: 
%                    w.........-> frequency 
%                    t.........-> time
%                    a.........-> amplitude
%                    bins .....-> number of the color map divisions
% Remarks..........:
%                    The colormap is ploted in 256 levels of gray

function [meanMHS,d_s,f] = DS(m_a_p,minFreq,maxFreq,hs_dt,T)

    [mhs,f] = MHSpec(m_a_p,minFreq,maxFreq,hs_dt);
    meanMHS = mhs/T;
    meanMHS = 10.^(10./meanMHS);
    Index = find(meanMHS==0);
    meanMHS(Index)=1e-20;
   
    [r,c]=size(m_a_p); 
   
    d_s = zeros(1,r);
   
    for k = 1:r,
        Index = find(m_a_p(k,:)==0);
        m_a_p(k,Index)=1e-20;
        d_s(1,k)=(trapz((1-(10.^(10./m_a_p(k,:)))./meanMHS(k)).^2).*hs_dt)./T;
    end;

