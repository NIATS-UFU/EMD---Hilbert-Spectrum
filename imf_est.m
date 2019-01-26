% Function name....: imf_est
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    imf_est estimates an intrinsic mode function (IMF). This function
%                    has the following characteristics: 
%                    (1) in the whole data set, the number of extrema and the number of 
%                        zero crossings must be either equal or differ at most by one; 
%                    (2) at any point, the mean value of the envelope defined by the local
%                        maxima and the envelope defined by the local minima is zero. It is 
%                        very difficult to obtain this requirement. Especially due to introduction
%                        of errors in the estimation of the upper and lower envelopes.
% Parameters.......: 
%                    x ......-> input time-series 
%                    IsIMF ..-> flag. If IsIMF>0 verify whether x is or not an IMF
%                    interpMethod ..-> interploation method: 0 = spline and ~0 pchip
%                    IMF_ID.........-> IMF identifier
% Return...........:
%                    i_m_f..-> intrinsic mode function
%                    is_imf.-> flag. If IsIMF>0 indicate that x is an IMF. Otherwise x is not an IMF!
% for compilation: mcc -x imf_est
function [i_m_f,is_imf]= imf_est(x,IsIMF,interpMethod,IMF_ID)

    %initializing variables
    i_m_f=[];
    i=0;
    stop_cond1 = 3;
    is_imf = -1;
    size_x = length(x);

    if(IsIMF>=0)
        %Before starting the sifting process it is interesting to check whether the input signal (x)
        %is an IMF or not.
        %Estimating the number of peaks and valleys of the signal
         [AveEnvelope,yu,yl,n_extrema]= meanEnv(x,interpMethod);
        %Estimating the number of zero crossings
         [z]= zero_crossing(x); 
        %Verifying the first condition
        if(abs(z-n_extrema)<=1 && sum(AveEnvelope)==0)
             i_m_f = x;
             disp('The input signal is an intrinsic mode function');
             fprintf(1,'Number of zero crossings: %d\n',z);
             fprintf(1,'Number of extreme: %d\n',n_extrema);
             fprintf(1,'Mean envelope: %f\n',sum(AveEnvelope)/length(AveEnvelope));
             disp('oi');
             is_imf = 1;
             return;
        end;%if
    end%if
    
    %Initiating sifting process
    while (stop_cond1 >= 2),

        i=i+1; %counter
        
        %Estimating the upper, lower and average envelope, and the number of extreme (peaks and valleys) 
        [AveEnvelope,yu,yl,n_extrema]= meanEnv(x,interpMethod);
        if (isempty(AveEnvelope)),
            return;
        end%if
        
        %Estimating the mean of AveEnvelope
        %This is important as a reference parameter that allows the user decide whether 
        %the estimated IMF will be used or not. This parameter is not being used as 
        %a stopping condition
        stop_cond2 = abs(sum(AveEnvelope)/size_x);
        fprintf(1,'IMF[%d]\n',IMF_ID);
        fprintf(1,'Average envelope: %f\n',stop_cond2);

        %Sifting the average envelope from the signal x
        h=x-AveEnvelope;
        %Checking whether h is or not an intrisic mode function
        %The difference between extrema and zeros crossings is less than or equal to 
        % two.
        [z]= zero_crossing(h); %number of zero crossings.
        stop_cond1 = abs(n_extrema-z);
        fprintf(1,'Zero-crossings - extreme =  %d\n',stop_cond1);
        %updating the input signal
        x=h;
    end;%while
    i_m_f = h;%updating the output variable
   