% Function name....: plotHS
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

function [mhs,f] = MHSpec(m_a_p,minFreq,maxFreq,hs_dt)

    [r,c]=size(m_a_p); 
    
    mhs = zeros(1,r);
    
    for k = 1:r,
         mhs(1,k)=trapz(m_a_p(k,:))*hs_dt;
    end

    f = 1:1:r;
    f = convScale(1,r,f,minFreq,maxFreq);
