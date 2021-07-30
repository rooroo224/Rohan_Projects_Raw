%% New modified code for residual plot
 
clear; clc; close all
clc;

%modifiable variables
n=7; l=n; nu1=2; nu2=1; gamma=2;

% Note here N-1 inside domain points are used along with 2 boundary points
N=2^n; Nc= N/2; h=1/N; H=1/Nc;
uh=zeros(N+1); ua=zeros(N+1); f=zeros(N+1);  
x=linspace(0,1,N+1); y=linspace(0,1,N+1);

% Calculate analytical solution ua and right hand side f
for i=1:N+1
    for j=1:N+1
         f(i,j)=8*pi^2*sin(2*pi*x(i))*sin(2*pi*y(j));
         ua(i, j)=sin(2*pi*x(i))*sin(2*pi*y(j));
    end
end

 global count 
 count=0;

 global a
 a=zeros(1,1);
 
 %Multigrid function
 tic
 [uout, count] = MG(l, uh, f, gamma, nu1, nu2, count);
 toc
 
%residual
tol = 10^(-10);m = 25;  
r_0 = Residual(uh,f);     %calculate initial residual
r0 = infnormcalc(r_0);
rel = zeros(m,1);   %store residual after every multigrid iteration
for i=1:m
    uh = MG(l, uh, f, gamma, nu1, nu2, count);
    r_m = Residual(uh,f);
    rm = infnormcalc(r_m);     % calculate infinite norm of vector r_m
    p = rm/r0;              
    rel(i) = p;
    if p < tol             % check relative residual
        break;
    end
end
 %Plotting
 figure (1)
 surf(uout)
 title('Numerical Solution Surface Plot with Multigrid Method');
 figure (2)
 surf(ua)
 title('Analytical Solution');


 figure(4)
 semilogy(1:i,rel(1:i))
 title('Çonvergence Plot');
 xlabel('Íterations');
 ylabel('log - InfNorm(rm/ro)');
 
 
 
 
 
 
 
 
 
 
 