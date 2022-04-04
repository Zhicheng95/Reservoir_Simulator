function [M_tp,grid_num]=matrix_tp(M,tag)
%%% please input the index matrix in the direction of N is down, same
%%% direction of i and E is right, the same direction of j.

% This function is used to transpose the matrix from the original direction
% into the target direction


% Pretreatment for input, transposing the matrix and flip it over.
[numj,numi]=size(M);

if tag==1
    [grid_num_ini]=numbering_final(M);
end

M_tp=zeros(numi,numj);
grid_num=zeros(numi,numj);

for i=1:numi
    for j=1:numj
        M_tp(i,j)=M(numj+1-j,numi+1-i);
        if tag==1
            grid_num(i,j)=grid_num_ini(numj+1-j,numi+1-i);
        end
    end
end

