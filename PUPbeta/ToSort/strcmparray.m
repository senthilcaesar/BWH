function [A,A1only,A2only,Both] = strcmparray(A1,A2)

for i=1:length(A1)
    for j=1:length(A2)
        A(i,j)=strcmp(A1(i),A2(j));
    end
end

Both = A1(sum(A,2)==1);
A1only = A1(sum(A,2)==0);
A2only = A2(sum(A)==0);
