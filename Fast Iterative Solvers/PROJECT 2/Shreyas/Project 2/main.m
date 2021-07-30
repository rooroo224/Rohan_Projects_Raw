clear

n = 7;
N = 2^n;                %total elements in each direction
h = 1/N;                %mesh size
u = zeros(N+1,N+1);     %store multigrid solution
f = zeros(N+1,N+1);     % rhs of poisson problem
u_true = zeros(N+1,N+1);
tol = 10^(-10);
m = 25;                 %approximate multigrid iterations

for i=1:1:N+1
    for j=1:1:N+1
        f(i,j) = 8*pi*pi*sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);       %store given RHS 'f'
        u_true(i,j) = sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);          % true solution for u
    end
end

s=cputime;
r_0 = Res(u,f);     %calculate initial residual
r0 = InfNorm(r_0);
rel = zeros(m,1);   %store residual after every multigrid iteration
for i=1:m
    u = MG( u, f, 2, 1, 1);
    r_m = Res(u,f);
    rm = InfNorm(r_m);     % calculate infinite norm of vector r_m
    p = rm/r0;              
    rel(i) = p;
    if p < tol             % check relative residual
        break;
    end
end
e=cputime-s;
semilogy(1:i,rel(1:i));
xlabel('Multigrid Iterations (m)');
ylabel('||r_m||_\infty /||r_0||_\infty');
title('Relative Residual v/s Multigrid Iterations');

% %plot solution u
% x = 0:h:1;
% y = 0:h:1;
% [x,y] = meshgrid(x,y);
% surf(x,y,u_true);
% zlim([-1 1]);
% title('Multigrid Solution, n = 4')
