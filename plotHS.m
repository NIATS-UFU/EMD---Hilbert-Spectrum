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

function plotHS(w,t,a,bins)

    sizeW = length(w);
    
    w_s = round(convScale(min(w),max(w),w,1,bins)); %color map row
    t_s = round(convScale(min(t),max(t),t,1,bins)); %color map column
    
    %creating color map
    m_a_p = min(a)*ones(bins);
    for i = 1:sizeW,
         m_a_p(w_s(i),t_s(i)) = a(i);
    end
    
    freq = 1:1:bins;
    freq = convScale(1,bins,freq,min(w),max(w));

    time = 1:1:bins;
    time = convScale(1,bins,time,min(t),max(t));
      
    pcolor(time,freq,m_a_p);
    shading('interp');
    colormap(gray(256));
 
