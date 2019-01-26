% Function name....: findClosest.m
% Date.............: October 15, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    Given two vectors, X1 and X2, this function generates a third
%                    vector, Z, obtained from a comparison between them.
%                    The main idea is to have the elements of X1 as
%                    reference and to search for the closest elements in
%                    X2. See the example bellow for a better understand:
%
%                    X1 = [1 2 3 4 5 6 7 8 9 10]; %this is the reference
%                    vector
%                    X2 = [3 6 9 15]; %this is the 
%                    Z = [3 3 3 ]

%
% Parameters.......: 
%                    x .....-> a scalar input 
% Return...........:
%                    s  -> outuput as defined above
% for compilation: mcc -x sng

function [Z]= findClosest(X1,X2)

Nsamples = length(X2);
for k=1:Nsamples,
    
    [amp,indx] = min(abs(X1 - X2(k)));
     Z(k) = X1(indx);
    
%      disp(k);
%     flag=1;
%     d=0; 
%     index = [];
%     while flag==1,
%         [index] = find(X1>=X2(k) & X1<=X2(k)+d);
%         d = d+1;
%         if isempty(index)~=1 | k==Nsamples | d>Nsamples,
%             flag= 0;
%         end;
%     end;
%     if isempty(index)~=1, Z(k) = X1(index(1)); end;
    
end;

   