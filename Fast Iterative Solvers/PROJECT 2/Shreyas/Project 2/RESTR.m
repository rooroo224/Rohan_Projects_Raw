function [u_H] = RESTR(u_h)

N = size(u_h,1)-1;
Nc = N/2;
h = 1/N;
H = 2*h;
u_H = zeros(Nc+1,Nc+1);

for i=2:1:Nc
    ii = 2*i-1;
    for j=2:1:Nc
        jj = 2*j-1;
        u_H(i,j) = (1/16)*(u_h(ii-1,jj-1) + 2*u_h(ii,jj-1) + u_h(ii+1,jj-1) + 2*u_h(ii-1,jj) + 4*u_h(ii,jj) + 2*u_h(ii+1,jj) + u_h(ii-1,jj+1) + 2*u_h(ii,jj+1) + u_h(ii+1,jj+1));
    end
end

end