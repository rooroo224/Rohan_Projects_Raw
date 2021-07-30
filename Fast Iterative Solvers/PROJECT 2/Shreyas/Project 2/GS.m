function [u] = GS(u,f,v)

N = size(u,1)-1;
h = 1/N;

for k=1:v
    for i=2:1:N
        for j=2:1:N
            u(i,j) = 0.25*(h*h*f(i,j) + u(i-1,j) + u(i,j-1) + u(i+1,j) + u(i,j+1));     % Gauss-Seidel step
        end
    end
end

end