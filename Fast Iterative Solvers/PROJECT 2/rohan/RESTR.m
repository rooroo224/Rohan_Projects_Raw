function [fnew]=RESTR(r,Nc)
fnew=zeros(Nc,Nc);
%Restriction Operator
for i=2:Nc-1
    ii=2*i;
    for j=2:Nc-1
        jj=2*j;
        fnew(i,j)=1/16 *(r(ii-1,jj-1)+ 2*r(ii,jj-1)+ r(ii+1,jj-1)+2*r(ii-1,jj)+4*r(ii,jj)+2*r(ii+1,jj)+r(ii-1,jj+1)+ 2*r(ii,jj+1)+ r(ii+1,jj+1));
    end
end  
end