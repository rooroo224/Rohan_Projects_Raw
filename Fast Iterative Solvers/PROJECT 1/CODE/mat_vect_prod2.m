function [Z]= mat_vect_prod2(values,colid,rowptr,rows,V)
%ref=0;k=0;
%mat-vect product for non-symmetric matrix using CRF
Z =zeros(rows,1);
for i=1:rows
    for k=rowptr(i):(rowptr(i+1)-1)
       Z(i,1)=Z(i,1)+values(k)*V(colid(k),1);
    end
end
% % Implementing Matrix vector product for Symmetric Matrices
% % Storing elements of Symmetric matric in CRS format
% if mat_type==1
%     for i=1:rows
%     for j=1:i
%         if (A(i,j)~=0)
%         ref=ref+1;
%         values(ref)= A(i,j);  
%         colid(ref)=j;
%         end
%     end
%     rowptr(i+1)=ref+1;
%     end
% %Symmetric Matrix Vector Product
%  Z=zeros(rows,1);
% for i=1:rows
%     for k=rowptr(i):(rowptr(i+1)-1)
%        Z(i,1)=Z(i,1)+ (values(k)*V(colid(k),1));
%     end
% end
%  for i=1:rows
%     for k=rowptr(i):(rowptr(i+1)-1)
%         Z(colid(k),1)=Z(colid(k),1)+values(k)*V(i,1);
%     end
%  end
%   p=zeros(1,rows);
% 
% for i=1:rows
%     Z(i,1)=Z(i,1)-values(rowptr(i+1)-1)*V(i,1);      
% end
% 
% end
% end





