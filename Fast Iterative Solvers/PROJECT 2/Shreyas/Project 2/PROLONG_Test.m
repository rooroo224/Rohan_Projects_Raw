
% Prolongation Operator - Work Package 2 (wp2)
clear all;

n = 4;
N = 2^n;
Nc = N/2;
h = 1/N;
H = 2*h;

u_h = zeros(N+1,N+1);
u_H = zeros(Nc+1,Nc+1);
u_true = zeros(N+1,N+1);      % true solution on Fine Grid

for i=1:1:Nc+1
    for j=1:1:Nc+1
        u_H(i,j) = sin(2*pi*(i-1)*H)*sin(2*pi*(j-1)*H);
    end
end

for i=1:1:N+1
    for j=1:1:N+1
        u_true(i,j) = sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);      % true solution on Fine Grid
    end
end

for i=2:1:Nc
    ii = 2*i - 1;
    for j=2:1:Nc
        jj = 2*j - 1;
        u_h(ii,jj) = u_h(ii,jj) + u_H(i,j); 
        u_h(ii-1,jj) = u_h(ii-1,jj) + 0.5*u_H(i,j);
        u_h(ii+1,jj) = u_h(ii+1,jj) + 0.5*u_H(i,j);
        u_h(ii,jj-1) = u_h(ii,jj-1) + 0.5*u_H(i,j);
        u_h(ii,jj+1) = u_h(ii,jj+1) + 0.5*u_H(i,j);
        u_h(ii-1,jj-1) = u_h(ii-1,jj-1) + 0.25*u_H(i,j);
        u_h(ii-1,jj+1) = u_h(ii-1,jj+1) + 0.25*u_H(i,j);
        u_h(ii+1,jj-1) = u_h(ii+1,jj-1) + 0.25*u_H(i,j);
        u_h(ii+1,jj+1) = u_h(ii+1,jj+1) + 0.25*u_H(i,j);
    end
end

max = 0;
for i=1:1:N+1
    for j=1:1:N+1
        if abs(u_h(i,j) - u_true(i,j)) > max
            max = abs(u_h(i,j) - u_true(i,j));
        end
    end
end

error_h = max;
        