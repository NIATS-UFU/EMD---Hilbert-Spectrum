% Function name....: random_sample
% Date.............: March 31, 2005
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    random sampling without replacement
% Parameters.......: 
%                    Population...........: vector with population samples
%                    Nsamples.............: number of required samples
% Return...........:
%                    y....................: vector with Nsamples samples sampled
%                                           from Population
function [y]=random_sample(Population,Nsamples)

n = length(Population);

if Nsamples>n,
    disp('The number of samples must be less than the size of the population.');
    y = [];
    return;
end;

p = randperm(n);
%permutation
Population = Population(p);
y = Population(1:Nsamples);