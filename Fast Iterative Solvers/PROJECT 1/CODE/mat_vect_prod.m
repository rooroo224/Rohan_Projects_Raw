function [Z]= mat_vect_prod(A,V)
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
mat_type=1; % mat_typr=0 for non-symmetric and else mat_type=1 for symmetric

%checking if imput matrix is symmetric or non symmetric
if rows~=cols
    mat_type=0;
elseif rows==cols
for i=1:rows  
    for j=1:cols
        if (A(i,j)~= A(j,i) || A(i,i)==0)
            mat_type=0;
        end
    end
end
end

% implementing mat-vect product for non-symmetric matrix
if mat_type == 0   
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
% Implementing Matrix vector product for Symmetric Matrices
% Storing elements of Symmetric matric in CRS format
if mat_type==1
    for i=1:rows
    for j=1:i
        if (A(i,j)~=0)
        ref=ref+1;
        values(ref)= A(i,j);  
        colid(ref)=j;
        end
    end
    rowptr(i+1)=ref+1;
    end
%Symmetric Matrix Vector Product
 Z=zeros(rows,1);
for i=1:rows
    for k=rowptr(i):(rowptr(i+1)-1)
       Z(i,1)=Z(i,1)+ (values(k)*V(colid(k),1));
    end
end
 for i=1:rows
    for k=rowptr(i):(rowptr(i+1)-1)
        Z(colid(k),1)=Z(colid(k),1)+values(k)*V(i,1);
    end
 end
  p=zeros(1,rows);

for i=1:rows
    Z(i,1)=Z(i,1)-values(rowptr(i+1)-1)*V(i,1);      
end

end
end





