% Function name....: Butfilter
% Date.............: November 12, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Low pass Butterworth filter
%                    
% Parameters.......: 
%                    signal .............: input signal
%                    fsr.................: sampling frequency
%                    fc..................: cut-off frequency
%                    n...................: order
% Return...........:
%                    outSignal...........: filtered signal
%                    ybf.................: FFT before filtering
%                    y...................: FFT after filtering


function [filterCoef,f,outSignal]= Butfilter (signal,fsr,fc,n,filterType)

filterCoef = [];
y = [];
ybf = [];
outSignal= [];

nrAmos = length(signal);

% h = waitbar(0,'Please wait...');
% waitbar(n/nrAmos);

FreqRes = fsr/nrAmos; %frequency resolution (Hz);

NsamplesFilter = (nrAmos-1)/2;

%Estimating filter coefficients 
% switch filterType,
%     case 1, %low-pass filter
        for i=1:(nrAmos-1)/2,
%             waitbar(f/(nrAmos/2));%updating waitbar
            % f and fsr (Hz)
            f(i) = i*FreqRes;
        	w = (2*pi*f(i))/fsr; %(rad/s)
        	filterCoef(i) = sqrt(1/(1+(((1/(tan(pi*fc/fsr)))*tan(w/2))^(2*n))));
         end;%for
%     case 2,
        

    
    
% end;


% 	case 1: //passa alta
% 			for(f=0;f<nrAmos;f++)
% 			{
% 				//Calculo da frequencia w (rad/s)
% 				//f e fsr (Hz)
% 				w = (2*pi*f)/fsr;
% 				filterCoef[f] = sqrt(1/(1+pow(((1/(1/tan(pi*fh/fsr)))*1/tan(w/2)),(2*N))));
% 			}//for f*/
% 			break;
% 	case 2: //rejeita faixa
% 			for(f=0;f<nrAmos;f++)
% 			{
% 				double aux1, aux2;
% 				//Calculo da frequencia w (rad/s)
% 				//f e fsr (Hz)
% 				w = (2*pi*f)/fsr;
% 				aux1 = sqrt(1/(1+pow(((1/(tan(pi*fl/fsr)))*tan(w/2)),(2*N))));
% 				aux2=sqrt(1/(1+pow(((1/(1/tan(pi*fh/fsr)))*1/tan(w/2)),(2*N))));
% 				if(aux2>aux1) filterCoef[f] = aux2;
% 				  else filterCoef[f] = aux1;
% 			}//for f



% % %Fast Fourier Transform of the input signal
y = fft(signal);
ybf=y;
  
k=1;
for i=2:(nrAmos-1)/2,
    y(i) = y(i)*filterCoef(k);
    k=k+1;
end;
  
k=1;
for i=nrAmos:-1:(nrAmos-1)/2+2,
    y(i) = y(i)*filterCoef(k);
    k = k+1;
end;
outSignal = real(ifft(y));

% close(h);%closing waiting bar



