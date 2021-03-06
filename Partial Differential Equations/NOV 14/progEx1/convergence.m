function error =convergence(u,x,f)
    y = linspace(0,1,1000);
    v=y;
    if f==1
        v =0.125-0.5*(y-0.5).^2;
    else
        for i = 1: size(y,2)
            if y(i)<0.5
                v(i)=y(i);
            else
                v(i) = 1-y(i);
            end
        end
    end

    u1 = interp1(x,u,y) ;
    error = sum((v - u1).^2);
end