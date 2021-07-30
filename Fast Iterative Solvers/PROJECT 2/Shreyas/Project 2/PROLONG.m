function [u_h] = PROLONG(u_H)

Nc = size(u_H,1)-1;
N = 2*Nc;
h = 1/N;
H = 2*h;
u_h = zeros(N+1,N+1);

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

end
