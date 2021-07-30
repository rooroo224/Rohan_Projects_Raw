function [u] = GS(u, f, nu)
[l,~]=size(u);
N=l-1;h=1/N;
%looping over internal variables
for k=1:nu  
    for i=2:N  
        for j=2:N  
            u(i,j)=1/4*(h^2*f(i,j)+u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1));
        end
    end
end
end

