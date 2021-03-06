clear;
clc;
L= 1; %length
N = 11:200:700;% no of internal nodes
E = zeros(size(N,2),1);
for i = 1:size(N,2)
    [x,h]= grid(N(i),L);
    f = 2;   %1 for f =1;
             %2 for f = dirac;
    [A,B]= assembly(x,h,f);
    u = linEqnsolver(A,B);
    e =convergence(u,x,f);
    E(i)= log(e);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
 plot(x,u);
 xlabel('x')
 ylabel('u')
 title('Solution of the PDE')
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Plot Convergence rate
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure;
plot(log(N),E)
xlabel('log(N)')
ylabel('log of Error')
title('CONVERGENCE RATE')