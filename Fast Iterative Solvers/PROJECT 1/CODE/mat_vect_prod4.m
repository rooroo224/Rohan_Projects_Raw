function [Z]= mat_vect_prod4(A,V)
[rows,cols] = size(A);
entries=0;
for i=1:rows
   for j=1:cols
       if A(i,j)~=0
       entries=entries+1;     
       end
   end
end
values=zeros(1,entries);colid=zeros(1,entries);
ref=0;k=0;rowptr(1)=1; colptr(1)=1; 
% implementing mat-vect product for non-symmetric matrix
% storing elements in CRS Format
for i=1:rows
    for j=1:cols
        if (A(i,j)~=0)
        ref=ref+1;
        values(ref)= A(i,j);  
        colid(ref)=j;
        end
    end
    rowptr(i+1)=ref+1;
end   
%mat-vect product for non-symmetric matrix using CRF
Z =zeros(rows,1);
for i=1:rows
    for k=rowptr(i):(rowptr(i+1)-1)
       Z(i,1)=Z(i,1)+values(k)*V(colid(k),1);
    end
end
end





