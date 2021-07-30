clear all;clc
%validating Restart GMRES
clear all; clc;
tic
 xo=zeros(1030,1);  %5357 1030
[A,rows,cols,entries,~,~,~] = mmread('orsirr_1.mtx');
xstar=ones(1030,1);v=A*xstar;  m=zeros(rows,cols); minv=zeros(rows,cols);
for i=1:rows
     for j=1:i
         m(i,j)= A(i,j);
     end
end
minv= inv(m);
%A=[0 5 0;0 2 3;1 0 7];  v=[10;-1;3]; xo=[0 ;1;0]; rows=3;cols=3;entries=5;
% storing elements in CRS Format
values=zeros(1,entries);colid=zeros(1,entries);rowptr(1)=1;ref=0;

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
[x,rho]=pre_gmres_gs(values,colid,rowptr,rows,v,xo,300,minv);
%[p]=mat_vect_prod2(values,colid,rowptr,rows,v)
%[x,rho,count,e]=restart_gmres3(values,colid,rowptr,rows,xo,v,100);
elapsedtime=toc
