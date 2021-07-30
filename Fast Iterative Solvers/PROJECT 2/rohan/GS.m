function [u]= GS(u,f,nu,N)
h=1/N;
    for k=1:nu
        for i=2:N-1
            for j=2:N-1
                u(i,j)=0.25*(h^2*f(i,j)+ u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1));
            end
        end
    end
end
