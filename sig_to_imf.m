% Function name....: sig_to_imf
% Date.............: December 13, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    sig_to_imf decomposes an input time-series x into its intrinsic mode functions.
% Parameters.......: 
%                    x     .....-> input time-series (row vector)
%                    mse   .....-> mean squared error 
%           interpMethod   .....-> if 0 then splines are used otherwise
%                                  pchip is employed
% Return...........:
%                    i_m_fs.....-> intrinsic mode functions(matrix: each
%                                  row is an IMF, where the first row is the first imf
%                                  and so on ...
%                    r..........-> this is the residue
% Remarks..........:
%                    The number of imfs generated is dependent on the mse. It is important 
%                    to observe that a very small value for the mse can produce functions
%                    that are not really imfs but a simple trend that will hardly contribute
%                    to the final spectrum of the original signal.
% Revision.........:
%                    1st: April/2003
% For compilation: mcc -x sig_to_imf (MatLab compiler should be available)
% Example........:
%                 [i_m_fs,r]= sig_to_imf(x,1e-5,1);
function [i_m_fs,r]= sig_to_imf(x,mse,interpMethod)

    %initializing variables
    t=length(x);
    imf_old=zeros(1,t);
    i_m_fs=[];%each line of this array is one imf
    r=x; %r:residue
    k=0; %iteration counter
    m_s_e = mse + 1;%mse: mean squared error
    
    %Displaying the chosen interpolation method
    if(interpMethod==0),
        disp('Interpolation method: splines');
    else
        disp('Interpolation method: pchip');
    end;
    
    %Before starting the sifting process is important to verify whether the input signal 
    % is or not an IMF
    disp('Verifying whether the input signal is an IMF or not');
    disp(' .......... ');
    disp(' ');  
    [i_m_f,is_imf]= imf_est(x,0,interpMethod,1);
    if is_imf>0,
        i_m_fs = i_m_f;
        return;
    else 
        disp(' The input signal is not an IMF and its first IMF has been estimated. ');
        imf_old = i_m_f;
        i_m_fs = i_m_f; %updating outuput vector
        k = k+1; %updating counter
    end;
    
    %Starting the sifting process
    while (m_s_e>mse)&(k<=9),

        k = k+1;%counter
        fprintf(1,'\nEstimating IMF: %d\n',k);        
        
        r=r-imf_old;
        [i_m_f,is_imf]= imf_est(r,-1,interpMethod,k);%estimating imf        
        if (isempty(i_m_f)),
            disp('A trend was detected');%there is not maximum or minimum local
            break;
%             return;
        end%if

        %updating the i_m_fs array
        i_m_fs = [i_m_fs; i_m_f];
        
        %Estimating the mean squared error between two consecutives imfs
        if (k>1),
            m_s_e=0;
            for i=1:t,
                m_s_e=m_s_e + (imf_old(i)-i_m_f(i))^2;
            end %for
        end%if
        m_s_e=m_s_e/t;
        fprintf(1,'\nmse: %d\n',m_s_e);        

        %updating old imf
        imf_old=i_m_f;
    end%while
    
    %Estimating the index of orthogonality
    auxImfs = [i_m_fs;r];
    numberOfIMFs = min(size(auxImfs));
    T = max(size(auxImfs));
    sqrX = x.^2;
    Index = find(sqrX==0);
    sqrX(Index)=1e-80;%very small value
    
    IO=0;%index of orthogonality
    for j=1:numberOfIMFs,
        sumk=0;
        for k=1:numberOfIMFs,
            sumk = sumk + sum(auxImfs(j,:).*auxImfs(k,:));
            sumk = sumk./sum(sqrX);
        end%for
        IO = IO + sumk;
    end%for
    fprintf(1,'\nOverall index of orthogonality: %f\n',IO);  

    
