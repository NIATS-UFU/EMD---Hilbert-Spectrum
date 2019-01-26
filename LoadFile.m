% Function name....: LoadFile
% Date.............: January 06, 2002
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    LoadFile loads a text file. The data stored in the file are float numbers 
%                    separated by spaces and with commas as a decimal symbol.
% Parameters.......: 
%                    filename.........-> file name with extension  
% Return...........:
%                    data.............-> vector containing the data loaded from the file.

function [data] = LoadFile(filename)

    data = [];

    fid = fopen(filename);%opening file
    
    if(fid==-1),%Is fid a valid handle?
        disp('Error: cannot open the file');
        return;
    end %if
    
    [A,count] =  fscanf(fid,'%c',inf);%reading file contents
    str=strrep(A,',','.'); %changing commas by periods
    data = sscanf(str,'%f')'; %converting string into float
    
    fclose(fid);%closing file
