% Function name....: getContourLevels
% Date.............: October 26, 2004
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This function gets contour levels from a matrix 'C' estimated 
%                    by either 'contour' or 'contour3'
% Parameters.......: 
%                    C.....................-> contours 
% Return...........: 
%                    contourLevels.........-> available contour levels

function [contourLevels] = getContourLevels(C)

counter = 1;
kk = 1;

contourLevels(kk) = C(1,counter); %first contour level
nElemts = C(2,counter); %number of points of the first patch
 disp( contourLevels(kk) );

while 1,
    
    counter = counter + nElemts + 1;%disp(counter);
    if(counter > size(C,2)), break; end; %Stop condition
    
    nElemts = C(2,counter);
    auxVar = C(1,counter);
    
    if(contourLevels(kk)~=auxVar)
        kk = kk + 1;
        contourLevels(kk) = C(1,counter);
        disp( contourLevels(kk) )
    end;
end;