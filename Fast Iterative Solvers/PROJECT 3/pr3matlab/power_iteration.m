close all; clear; clc;
N_Iter=1000; tol = 10^-8;  %Number of Power Iterations and Tolerance

% Pure Power Iteration
%Reading the Matrix in .mtx format and converting into CSR

[A,rows,cols,entries,~,~,~] = mmread('nos6.mtx');
%[A,rows,cols,entries,~,~,~] = mmread('s3rmt3m3.mtx');

tic
values=zeros(1,entries);colid=zeros(1,entries);
ref=0;rowptr(1)=1;
for i=1:rows
    for j=1:i
       if (A(i,j)~=0)
          ref=ref+1;
          values(ref)= A(i,j);  
          colid(ref)=j;
       end
    end
    rowptr(i+1) = ref+1;  
 end
  
%Implementing the Power Iteration Method
q = ones(rows,1)*1/sqrt(rows); tq = zeros(1,rows); lamda_old=0; diff = zeros(N_Iter,1);

q = mat_vect_prod3(values,colid,rowptr,rows,q);temp1 = q/norm(q,2); last=N_Iter;
for i=1:N_Iter
    q = temp1;
    temp = mat_vect_prod3(values,colid,rowptr,rows,temp1);
    temp1 = temp/norm(temp,2);
    %transpose
    for p = 1:rows
        tq(1,p) = q(p,1);
    end
    
    % q^T*A*q
    lamda = tq*temp;
    diff(i) = abs(lamda-lamda_old);
     
    if (diff(i)) < tol
        last = i;
        break;
    end
    lamda_old = lamda;
end
 toc
 formatSpec1 = 'The Approximate largest Eigan Value with Pure Power Iteration is %f\n';
 fprintf(formatSpec1,lamda)
 
 error = diff(last);
 
 formatSpec2 = 'The Error with Pure Power Iteration is %.15f\n';
 fprintf(formatSpec2,error)
 
 figure(1)
 semilogy(diff);
 title('Convergence Plot for Pure Power Iteration');
 xlabel('Iteration index'); ylabel('|\lambda^{k} - \lambda^{k-1}|');
 

