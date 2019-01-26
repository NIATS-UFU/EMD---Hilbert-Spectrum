% Function name....: threshold
% Date.............: April 23, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Automatic threshold estimation of EMG signals
% Parameters.......: 
%                    y .........-> input time-series
%                    Nbins .....-> number of bins for histogram estimation. For real signals it is related to the number of 
%                                  bits of the acquisition board (e.g. 2^(number_of_bits))
% Return...........:
%                    n...........-> fequency (number of occurrences)
%                    xout........-> time-series amplitude
%                    m...........-> mean of the distribution
%                    standardDev.-> standard deviation of the distribution

function [n,xout,m,standardDev,th,yenv]= threshold(y,Nbins)

%Estimating the complex envelope of the time series based on the hilbert transform
yenv = abs(hilbert(y));

%Estimating the histogram
[n,xout] = hist(yenv,Nbins); %number of bins = 2^9

% If the probability density function of the original signal, y, may be approximated by a Gaussian distribution
% then the peak of this distribution will be the mean (that is equal to the mode in this particular case) and 
% the standard deviation of this distribution can be obtained from a property of Gaussian distributions: 
% At x = s, where s is the standard deviation of the Gaussian distribution, the Gaussian falls 1/(sqrt(exp(1)) = 0.6065
% of its maximum value (aproximatelly 61%).

%Estimating the peak of the distribution (mean, m)
[m,Index_m] = max(n);

%Estimating the standard deviation
[index] = find(n(Index_m:end)<=0.61*m);
standardDev = n(index(1));

%Estimating threshold
th = 5*xout(index(1)); %based on the Chebychev's theorem