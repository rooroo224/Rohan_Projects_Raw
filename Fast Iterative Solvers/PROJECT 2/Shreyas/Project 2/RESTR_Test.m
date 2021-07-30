
% Restriction Operator - Work Package 2 (wp2)

clear all;

n = 4;
N = 2^n;
Nc = N/2;
h = 1/N;
H = 2*h;
u_h = zeros(N+1,N+1);
u_H = zeros(Nc+1,Nc+1);
u_true = zeros(Nc+1,Nc+1);      % true solution on Coarse Grid

for i=1:1:N+1
    for j=1:1:N+1
        u_h(i,j) = sin(2*pi*(i-1)*h)*sin(2*pi*(j-1)*h);
    end
end

for i=2:1:Nc
    ii = 2*i-1;
    for j=2:1:Nc
        jj = 2*j;
        u_H(i,j) = (1/16)*(u_h(ii-1,jj-1) + 2*u_h(ii,jj-1) + u_h(ii+1,jj-1) + 2*u_h(ii-1,jj) + 4*u_h(ii,jj) + 2*u_h(ii+1,jj) + u_h(ii-1,jj+1) + 2*u_h(ii,jj+1) + u_h(ii+1,jj+1));
        u_true(i,j) = sin(2*pi*(i-1)*H)*sin(2*pi*(j-1)*H);     %store true solution
    end
end

max=0;
for i=1:1:Nc+1
    for j=1:1:Nc+1
        if abs(u_H(i,j) - u_true(i,j)) > max
            max = abs(u_H(i,j) - u_true(i,j));
        end
    end
end
error_2h = max;        
        
        