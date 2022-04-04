function [grid_num]=numbering_final(index)
% This function is used to numbering the reservoir grid based on the rule
% of Standard ordering by columns

% This numbering is only for the final, left direction is North, and i direction,
% and up is East, and j direction.
% x direction is west direction and y direction is 
[numi,numj]=size(index);
grid_num=zeros(numi,numj);
count=1;
for i=1:numi
    for j=1:numj
        if index(numi+1-i,numj+1-j)~=0
            grid_num(numi+1-i,numj+1-j)=count;count=count+1;
        end
        
    end
end

            
