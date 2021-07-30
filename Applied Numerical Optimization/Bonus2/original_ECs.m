function [c,ceq] = original_ECs(n, m, Q, A, b, x)
% Define equality constraints for original problem to use fmincon solver.
c = [];
ceq = zeros(m,1);
for k = 1:m
    Qi = Q(n*(k-1)+1: n*k, :);
    ceq(k) = x'* Qi * x + A(k,:) * x - b(k); 
end

end

