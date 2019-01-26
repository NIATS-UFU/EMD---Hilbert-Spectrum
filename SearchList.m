% Function name....: SearchList
% Date.............: July 29, 2003
% Author...........: Adriano de Oliveira Andrade
% Description......:
%                    This function looks for elements of VecList2 in VecList1 with tolerance
%                    especified by Tolerance
% Parameters.......: 
%                    VecList1 ........-> Lookup table where values of VecList2 will be sought 
%                    VecList2 ........-> Reference amplitudes for the search
%                    Tolerance........-> Query resolution (it depends upon the units of VecList1 and VecList2)
% Return...........:
%                    QueryResult .... -> It returns indexes of VecList1 which are assigned to VecList2

function [QueryResult]= SearchList(VecList1,VecList2,Tolerance)

    %Initializing variables 
    sizeVecList2 = length(VecList2);
    QueryResult = -1*ones(1,sizeVecList2);% -1 means there is not assigned index 
    
    for i=1:sizeVecList2,
        
        [index] = find(VecList1>=VecList2(i)-Tolerance & VecList1<=VecList2(i)+Tolerance);
        if (index~=[]),
            MinIndex = min(index);
            QueryResult(i) = MinIndex;
        end %if
    end%for

    %Query return
    INDEX = find(QueryResult>0);
    AssignedIndex = QueryResult;
    QueryResult = VecList1(QueryResult(INDEX));
    
    k=1;
    for i=1:length(AssignedIndex),
        if(AssignedIndex(i)>0),
            AssignedIndex(i) = QueryResult(k);
            k = k+1;
        end%if
    end%for
    
    QueryResult = AssignedIndex;
        
    
    
    
    
    