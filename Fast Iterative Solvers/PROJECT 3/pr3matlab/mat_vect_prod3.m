function [Z]= mat_vect_prod3(values,colid,rowptr,rows,V)
% Implementing Matrix vector product for Symmetric Matrices
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





