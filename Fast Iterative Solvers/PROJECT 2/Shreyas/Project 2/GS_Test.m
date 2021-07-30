
% Gauss Seidel Relaxation - Work Package 1 (wp1)

clear all;

N = 100;
iter = 10000;
h = 1/N;
f = zeros(N+1,N+1);
u = zeros(N+1,N+1);
u_true = zeros(N+1,N+1);
tol = 10^(-10);

for i=1:1:N+1
    for j=1:1:N+1
        f(i,j) = 8*pi*pi*sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);       %store given RHS 'f'
        u_true(i,j) = sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);
    end
end

for k=1:iter
    w = u;
    for i=2:1:N
        for j=2:1:N
            u(i,j) = 0.25*(h*h*f(i,j) + u(i-1,j) + u(i,j-1) + u(i+1,j) + u(i,j+1));     % Gauss-Seidel step
        end
    end
    p = max(u-w,[],'all');
    if p < tol
        break;
    end
end

max = 0;
for i=1:1:N+1
    for j=1:1:N+1
        if abs(u(i,j) - u_true(i,j)) > max
            max = abs(u(i,j) - u_true(i,j));
        end
    end
end

error_max = max;

        