% Backward heat equation in 1D
% Bonus problem 3: iterative regularization + parameter choice

close all;
clear all;

L = pi;       % length of the rod
T = 5;        % final time
sigma = 1e-2; % magnitude of error
kappa = 1e-3; % heat conduction coefficient
n = 100;
N = n+1;      % number of intervals
h = L/N;      % grid size

% create C
C=zeros(n,n);
% first and last row
C(1,1:2)=[2 -1];
C(n,n-1:n)= [-1 2];
% inner rows
for i=2:n-1
    C(i,i-1:i+1)=[-1 2 -1];
end
% and, finally, multiply with the right factor
C=C*(kappa/h^2);

%% backward heat conduction
x=linspace(h,L-h,n)';       % generate linearly spaced vector
% u_0= sin(x)+0.2*sin(8*x); % initial temperature
u_0= exp(-(x-L/2).^2/0.02); % initial temperature
A= expm(-T*C);
u_f= A*u_0;                 % final temperature
load ('C:\Users\Rohan Krishna\Documents\Simulation Science\SEM 2\MBEM\BONUS\bonus3\error.mat');             % error (random noise)
u_f_meas= u_f + error;      % measurement of final temperature
u_0_est = expm(T*C)*u_f_meas;

sprintf('1d backward heat equation with n=%i and noise level sigma=%g\n==============================================================', n, sigma)
  
err= norm( u_0 - u_0_est)

figure
subplot(1,2,1);
plot(x, u_f, x, u_f_meas, x, A*u_0_est);
legend('u_f', 'u_f measured', 'u_f est. (data fit)');
subplot(1,2,2);
plot(x, u_0, x, u_0_est, 'r');
legend('u_0', 'estimation for u_0');

%% Landwever

% Landweber iterations
max_iter_LW= 250;

iter_LW= zeros(n,max_iter_LW); % store all iterates

% my code...
discrep_LW= zeros(1,max_iter_LW);
sol_norm_LW= zeros(1,max_iter_LW);
k_LW = [1:max_iter_LW];
beta_LW= 1;

ATA= A'*A;
ATym= A'*u_f_meas;
discrep_LW(1)= norm(0 - u_f_meas);
for i = 2 : max_iter_LW
    iter_LW(:,i) = iter_LW(:,i-1) - beta_LW * ATA * iter_LW(:,i-1) + beta_LW * ATym;
    discrep_LW(i)= norm(A*iter_LW(:,i) - u_f_meas);
    sol_norm_LW(i)= norm(iter_LW(:,i));
   
end

% plot figure
figure
data_error= 0.1;
subplot(1,2,1);
plot( k_LW, data_error*ones(size(k_LW)), k_LW, discrep_LW);
legend( '|| \delta y ||_2', '|| A \xi_{\alpha} - y_{meas} ||_2');
title( 'LW: Discrepancy');
subplot(1,2,2);
loglog( discrep_LW, sol_norm_LW, 'o-');
xlabel('||A \xi_{\alpha} - y_{meas}||_2');
ylabel('||\xi_{\alpha}||_2');
title( 'LW: L-curve');
for i=1:10:length(k_LW) % add text
    text(discrep_LW(i), sol_norm_LW(i), ['k = ' num2str(k_LW(i))]);
end

%pick k = 11 for discrepency principle and plot
figure
subplot(1,2,1);
plot(x, u_f, x, u_f_meas, x, A*iter_LW(:,11));
legend('u_f', 'u_f measured', 'u_f est. (data fit)');
title('Discrepancy for LW, k=11: data')
subplot(1,2,2);
plot(x, u_0, x, iter_LW(:,11), 'r');
legend('u_0', 'estimation for u_0');
title('Discrepancy for LW, k=11: exact and estimated parameters')

%pick k = 96 for L-curve and plot
figure
subplot(1,2,1);
plot(x, u_f, x, u_f_meas, x, A*iter_LW(:,96));
legend('u_f', 'u_f measured', 'u_f est. (data fit)');
title('L-curve for LW, k=96: data')
subplot(1,2,2);
plot(x, u_0, x, iter_LW(:,96), 'r');
legend('u_0', 'estimation for u_0');
title('L-curve for LW, k=96: exact and estimated parameters')

%% CGNE

% CGNE iteration
max_iter_CG= n;

iter_CG= zeros(n,max_iter_CG); % store all iterates

% my code...
discrep_CG= zeros(1,max_iter_CG);
sol_norm_CG= zeros(1,max_iter_CG);

k_CG = [1:max_iter_CG];
r = u_f_meas - 0; %r_0
d = A' * r;       %d_0
discrep_CG(1)= norm(0 - u_f_meas); %discrep_0
for i = 2 : max_iter_CG
    beta_CG= (norm(A'*r) * norm(A'*r)) / (norm(A*d) * norm(A*d));   
    iter_CG(:,i)= iter_CG(:,i-1) + beta_CG * d;
    rq= r; % save r_k-1
    r= r - beta_CG * A * d; % update r_k
    gamma= (norm(A'*r) * norm(A'*r)) / (norm(A'*rq) * norm(A'*rq)); %update gamma_k
    d= A' * r + gamma * d; %update d
    
    discrep_CG(i)= norm(A*iter_CG(:,i) - u_f_meas);
    sol_norm_CG(i)= norm(iter_CG(:,i));
   
end

% plot figure
figure
data_error= 0.1;
subplot(1,2,1);
plot( k_CG, data_error*ones(size(k_CG)), k_CG, discrep_CG);
legend( '|| \delta y ||_2', '|| A \xi_{\alpha} - y_{meas} ||_2');
title( 'CGNE: Discrepancy');
subplot(1,2,2);
loglog( discrep_CG, sol_norm_CG, 'o-');
xlabel('||A \xi_{\alpha} - y_{meas}||_2');
ylabel('||\xi_{\alpha}||_2');
for i=1:10:length(k_CG) % add text
    text(discrep_CG(i), sol_norm_CG(i), ['k = ' num2str(k_CG(i))]);
end
title( 'CGNE: L-curve');

%pick k = 5 for discrepency principle
figure
subplot(1,2,1);
plot(x, u_f, x, u_f_meas, x, A*iter_CG(:,5));
legend('u_f', 'u_f measured', 'u_f est. (data fit)');
title('Discrepancy for CGNE, k=5: data')
subplot(1,2,2);
plot(x, u_0, x, iter_CG(:,5), 'r');
legend('u_0', 'estimation for u_0');
title('Discrepancy for CGNE, k=5: exact and estimated parameters')

%pik k =21 for L-curve
figure
subplot(1,2,1);
plot(x, u_f, x, u_f_meas, x, A*iter_CG(:,21));
legend('u_f', 'u_f measured', 'u_f est. (data fit)');
title('L-curve for CGNE, k=21: data')
subplot(1,2,2);
plot(x, u_0, x, iter_CG(:,21), 'r');
legend('u_0', 'estimation for u_0');
title('L-curve for CGNE, k=21: exact and estimated parameters')

