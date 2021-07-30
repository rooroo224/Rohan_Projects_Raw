clear all; clc;
tic
%converting coo filr to comressed row format
[A,rows,cols,entries,~,~,~] = mmread('s3rmt3m3.mtx');
values=zeros(1,entries);colid=zeros(1,entries);
ref=0;k=0;rowptr(1)=1; 
%v=ones(5357,1);
% Storing elements of Symmetric matric in CRS format
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
 % [p]=mat_vect_prod3(values,colid,rowptr,rows,v)    
  
   xo=zeros(5357,1);  %5357 1030
   xstar=ones(5357,1);
   v=A*xstar;
   [d,e,count]=cg2(values,colid,rowptr,rows,v,xo);
   elapsedtime=toc
   %ploting
   for i=1:count
   cc= A*e(:,i);
   dd= dot(cc,e(:,i));
   nor(i,1)= sqrt(dd);
   end
    figure(2)
    x=1:count
    semilogy(x,nor)
    xlabel('Number of Iterations');
    ylabel('A norm (semilog)');
    title('Congugate Gradient-||e||_A');
   
  
   