%Sijie Luo            -403598
%Rohan Krishna Balaji -403596

%clear all;
clc;

% Reading Measurements of uc
uc= csvread('uc.csv');  

%initializing time, initial condition and parameters
t_tilda= 0:1:200; x0= [0;1];
c=2; L=12; R1=0.5; R2=4;

%Expressing the Model in Continuous LTI System
A=[-1/(c*R2) 1/c; -1/L R1/L]; B=[0; 1/L]; C=[1 0]; D=[0 0];

%Defining intercept function handle
uc_interp=@(t)interp1(t_tilda, uc, t);

% Choosing negetive poles
% Note different poles can be experimented to find optimal value,
% In general for "Larger Poles" for the estimated output fits well with the measured output 
P1=[-0.1 -0.2]; 
%For this "smaller p" case estimated current doesnot show high frequent oscillations but voltage doesnot fit well 

P2=[-1 -2];
%For this "medium p" case estimated current shows oscillations but voltage almost fits with measured

P3=[-5 -10];
%For this "medium p" case estimated current shows oscillations but voltage almost fits with measured

%Using pole placement function from matlab 
K1= place(A',C',P1);
K2= place(A',C',P2);
K3= place(A',C',P3);

%Defining Input Function
u=@(t) (5*(1-exp(-0.5*t)));

%Defing RHS for the Luenberger observer 
rhs_L1= @(t,x1) ((A-K1'*C)*x1+ B*u(t)+K1'*uc_interp(t));
rhs_L2= @(t,x2) ((A-K2'*C)*x2+ B*u(t)+K2'*uc_interp(t));
rhs_L3= @(t,x3) ((A-K3'*C)*x3+ B*u(t)+K3'*uc_interp(t));
%Solving the Differential equation
[~,x1]=ode45(rhs_L1,t_tilda,x0);
[~,x2]=ode45(rhs_L2,t_tilda,x0);
[~,x3]=ode45(rhs_L3,t_tilda,x0);

%Finding the exact state variables
rhs= @(t,x_hat) (A*x_hat +B*u(t));
[~, x_hat]=ode45(rhs,t_tilda,x0);

%Plotting
figure(1)
subplot(3,2,1);
plot(t_tilda,x1(:,1),'k-');
hold on;
plot(t_tilda,x_hat(:,1),'b-')
hold on;
plot(t_tilda,uc, 'r-');
legend('Uc estimates','Uc exact','measured')
title('Voltage (Uc) State Variable Small P');
xlabel('time'); ylabel('voltage');

subplot(3,2,3);
plot(t_tilda,x2(:,1),'k-');
hold on;
plot(t_tilda,x_hat(:,1),'b-')
hold on;
plot(t_tilda,uc, 'r-');
legend('Uc estimates','Uc exact','measured')
title('Voltage (Uc) State Variable Medium P');
xlabel('time'); ylabel('voltage');

subplot(3,2,5);
plot(t_tilda,x3(:,1),'k-');
hold on;
plot(t_tilda,x_hat(:,1),'b-')
hold on;
plot(t_tilda,uc, 'r-');
legend('Uc estimates','Uc exact','measured')
title('Voltage (Uc) State Variable Large P');
xlabel('time'); ylabel('voltage');

subplot(3,2,2);
plot(t_tilda,x1(:,2),'g-');
hold on;
plot(t_tilda,x_hat(:,2),'r-')
legend('I estimates','I exact');
title('Current(I) State Variable Small P');
xlabel('time'); ylabel('current');

subplot(3,2,4);
plot(t_tilda,x2(:,2),'g-');
hold on;
plot(t_tilda,x_hat(:,2),'r-')
legend('I estimates','I exact');
title('Current(I) State Variable Medium P');
xlabel('time'); ylabel('current');

subplot(3,2,6);
plot(t_tilda,x3(:,2),'g-');
hold on;
plot(t_tilda,x_hat(:,2),'r-')
legend('I estimates','I exact');
title('Current(I) State Variable Large P');
xlabel('time'); ylabel('current');






