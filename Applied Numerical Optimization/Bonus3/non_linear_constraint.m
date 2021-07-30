function [c,ceq] = non_linear_constraint(y,m,nx,nu)
        n = nx+nu;
        tf = 5; to = 0;
        h = ((tf-to)/(m));
        ceq = zeros(3*m,1);
        c = [];
        % initial conditons
        x1 = 0 ; x2 = 1; x3 = 0; u0 = 0;
       
        %non linear constraints 
        ceq(1,1) = y(1,1) - x1 - h*(((1 - y(2,1)^2)*y(1,1)) - y(2,1) + y(4,1));
        ceq(2,1) = y(2,1) - x2 - h*y(1,1);
        ceq(3,1) = y(3,1) - x3 - h*( y(1,1)^2 + y(2,1)^2 + y(4,1)^2 );
        
        % Implicit Eular Discretization
        count = nx;
        index = 1;
        for i = 2:m
            index = index + n;
            count = count + 1;
            ceq(count,1) = y(index,1) - y(index-n,1) - h*(((1-y(index+1,1)^2)*y(index,1)) - y(index+1,1) + y(index+3,1)); 
            count = count + 1;
            ceq(count,1) = y(index+1,1) - y(index-(n-1),1) - h*(y(index,1));
            count = count + 1;
            ceq(count,1) = y(index+2,1) - y(index-(n-2),1) - h*( y(index,1)^2 + y(index+1,1)^2 + y(index+3,1)^2 );
        end
end