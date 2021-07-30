function [obj_val] = definitions(m)
    nx = 3 ; nu = 1;
    n = nx+nu;

    % Lower Bound
    lb = -Inf*ones(n*m,1);

    for i = 1:n:n*m
        lb(i,1) = -0.4;
        lb(i+3,1) = -0.3;
    end

    % Upper Bound
    ub = Inf*ones(n*m,1);
    for i = 1:n:n*m
        ub(i+3,1) = 1;
    end

    % Objetive Function
    func = @(y) y((n*m)-1,1);
    
    % Equality and Inequality Constraints
    Aeq = [];
    beq = [];
    A = [];
    b = [];

    % Initializarion
    yo = ones(n*m,1);
    
    % Non linear constraint function definition
    nlc = @(y)non_linear_constraint(y,m,nx,nu);

    % Non Linear Optimization
    options = optimoptions('fmincon','MaxFunctionEvaluations',100000);
    y = fmincon(func,yo,A,b,Aeq,beq,lb,ub,nlc,options);

    % Objective Value
    obj_val = func(y);
end