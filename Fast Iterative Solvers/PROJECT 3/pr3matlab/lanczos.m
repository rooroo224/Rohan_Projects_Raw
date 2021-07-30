 %% Plot for Error v/s Time (Hard Coded from data collected)
 % If Pure iteration is needed uncomment the third section
 % If convergence is needed uncomment the fig 1 and 2 respectively
 
 clear; close all; clc; 
 tlan = [0.002316 ,0.005030 ,0.006669 ,0.025214];
 tpure = [0.060369, 0.393546, 1.069396 ,3.228686];
 
 e1 = [0.009672274087279,9.944619523594156*10^-05,9.970844985218719*10^-07,9.82254277914763*10^-11];
 e2 = [0.009961754039978,9.953436347132083*10^-05,9.993927960749716*10^-07,9.822542779147625*10^-11];
 
 subplot(2,2,1);
 semilogy(tlan,e1,'-o r');
 title('Error vs Time for Lancoz Method');
 xlabel('Time (in sec)'); ylabel('Error (Semilog)');
 
 subplot(2,2,2);
 semilogy(tpure,e2, '-o b');
 title('Error vs Time for Pure Power Iteration Method');
 xlabel('Time (in sec)'); ylabel('Error (Semilog)');

%%
% Implementing the Lancozs Method
% Run the first section if only Lanczos method is needed, if full run is
% performed both the Lancozs and Pure power methdods are executed and plotted

N_Iter=5000; %Number of Power Iterations
m = 100;   %Number of Krylov Spaces  
tol = 10^-10; %Tolerance

fprintf('Lanczos Method : ');
%Reading the Matrix in .mtx format and converting into CSR
[A,rows,cols,entries,~,~,~] = mmread('s3rmt3m3.mtx');

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

% Lanczos Method

v = zeros(rows,m+2); beta_temp = zeros(m+1,1); alpha = zeros(m,1);

for i = 1:rows
    v(i,2) = 1/sqrt(rows); %initialization
end

for j = 2:m+1
    
    w = mat_vect_prod3(values,colid,rowptr,rows,v(:,j)) - beta_temp(j-1,1)*v(:,j-1);
    alpha(j-1,1) = v(:,j)' * w;
    w = w - (alpha(j-1,1) * v(:,j));
    beta_temp(j,1) = norm(w,2);
    v(:,j+1) = w/beta_temp(j,1);
    
end

beta =  zeros(m-1,1);
for i = 1:size(beta_temp,1)-2
    beta(i,1) = beta_temp(i+1,1);
end

values = zeros(m+(m-1)*2,1); colid = zeros(m+(m-1)*2,1);
rowptr = zeros(m+1,1);
adder=0;
for i = 1:m-1
        values(adder+1) = alpha(i);
        colid(adder+1) = i;
        values(adder+2) = beta(i);
        colid(adder+2) = i+1;
        values(adder+3) = beta(i);
        colid(adder+3) = i;
        adder = adder+3;
end

values(m+(m-1)*2) = alpha(m);
colid(m+(m-1)*2) = m;
  
adder = 3;
rowptr(1,1)=1;
for i=2:m
   rowptr(i) = adder;
   adder = adder+3;
end
rowptr(m+1,1) = rowptr(m,1)+2;

%Power Iteration with reduced size
rows=m;
q = ones(rows,1)*1/sqrt(rows); tq = zeros(1,rows); lamda_old=0; diff = zeros(N_Iter,1); last = N_Iter;

tic
q = mat_vect_prod2(values,colid,rowptr,rows,q);temp1 = q/norm(q,2);
for i=1:N_Iter
    q = temp1;
    temp = mat_vect_prod2(values,colid,rowptr,rows,temp1);
    temp1 = temp/norm(temp,2);
    
    %transpose
    for p = 1:rows
        tq(1,p) =  q(p,1);
    end
    
    % q^T*A*q
    lamda = tq*temp;
    diff(i) = abs(lamda-lamda_old);
     
    if (diff(i)/diff(1)) < tol
        last = i;
        break;
    end
    lamda_old = lamda;
end
 toc
 formatSpec = 'The Approximate largest Eigan Value with Lanczos Iteration is %f\n';
 fprintf(formatSpec,lamda)
 
 error = diff(last);
 
 formatSpec2 = 'The Error with Lanzcos Iteration is %.15f\n';
 fprintf(formatSpec2,error)
 
%  figure(1)
%  semilogy(diff);
%  title('Convergence Plot for Lanczos Power Iteration');
%  xlabel('Iteration index'); ylabel('|\lambda^{k} - \lambda^{k-1}|');
%  
  %% Pure Power Iteration Method for full matrix
% %close all; clear; clc 
% %N_Iter=5000;          %uncomment only when running this section alone
% %tol = 10^-8;              %tolerance for pure power iteration 
% fprintf('Full Matrix Pure Power Iteration : ');
% 
% % Pure Power Iteration
% %Reading the Matrix in .mtx format and converting into CSR
% [A,rows,cols,entries,~,~,~] = mmread('s3rmt3m3.mtx');
% 
% 
% values=zeros(1,entries);colid=zeros(1,entries);
% ref=0;rowptr(1)=1;
% for i=1:rows
%     for j=1:i
%        if (A(i,j)~=0)
%           ref=ref+1;
%           values(ref)= A(i,j);  
%           colid(ref)=j;
%        end
%     end
%     rowptr(i+1) = ref+1;  
%  end
%   
% %Implementing the Power Iteration Method
% tic
% q = ones(rows,1)*1/sqrt(rows); tq = zeros(1,rows); lamda_old=0; diff = zeros(N_Iter,1);
% 
% q = mat_vect_prod3(values,colid,rowptr,rows,q);temp1 = q/norm(q,2); last = N_Iter;
% for i=1:N_Iter
%     q = temp1;
%     temp = mat_vect_prod3(values,colid,rowptr,rows,temp1);
%     temp1 = temp/norm(temp,2);
%     
%     %transpose
%     for p = 1:rows
%         tq(1,p) =  q(p,1);
%     end
%     
%     % q^T*A*q
%     lamda = tq*temp;
%     diff(i) = abs(lamda-lamda_old);
%      
%     if (diff(i)/diff(1)) < tol
%         last = i;
%         break;
%     end
%     lamda_old = lamda;
% end
%  toc
%  formatSpec = 'The Approximate largest Eigan Value with Pure Power Iteration is %f\n';
%  fprintf(formatSpec,lamda)
%  
%  error2 = diff(last);
%  
%  formatSpec2 = 'The Error with Pure Power Iteration is %.15f\n';
%  fprintf(formatSpec2,error2)
%  
% %  figure(2)
% %  semilogy(diff);
% %  title('Convergence Plot for Pure Power Iteration');
% %  xlabel('Iteration index'); ylabel('|\lambda^{k} - \lambda^{k-1}|');
 


 
 
 
 
 
 
 
 
 