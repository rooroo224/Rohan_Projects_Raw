function [r] = Res(u,f)

N = size(u,1)-1;
h =1/N;
r = zeros(N+1, N+1);

for i=2:1:N
    for j=2:1:N
        r(i,j) = f(i,j) + (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4*u(i,j))/(h*h);
    end
end

end