function [f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)
%Sijie Luo 403598
%Rohan Krishna 403596
%   This function is used to construct a convex relaxation of bilinear
%   terms by McCormick Envelops.
%   Returns the lower bound f_lb and an upper bound f_ub on the optimal
%   value of the optimization problem.

% 1. Check the inputs for correct length. Qi should be symmetric.
classes = {'numeric'};
validateattributes(c,classes,{'size',[n, 1]});     % check for dimensions
validateattributes(H,classes,{'size',[n, n]});
validateattributes(Q,classes,{'size',[m*n, n]});
validateattributes(A,classes,{'size',[m, n]});
validateattributes(b,classes,{'size',[m, 1]});
validateattributes(lb,classes,{'size',[n, 1]});
validateattributes(ub,classes,{'size',[n, 1]});
for k = 1:m
    if ~issymmetric(Q(n*(k-1)+1: n*k, :))
        error('Qi is not symmetric');
    end
end
disp('All inputs are fine.')

% 2. Analyze Q for how many auxiliary variables w need to be generated.
idx_w = []; % Store index of w
for k = 1:m
    Qi = Q(n*(k-1)+1: n*k, :);
    idx_w = [find(triu(Qi, 1)); idx_w]; 
end
uidx_w = unique(idx_w); % Unique index of w
n_w = length(uidx_w);   % The number of auxiliary variables w

% 3. Generate matrix B_ic for implementing the ICs.
% Variables becomes x = [x1, x2, ..., xn, w1, w2, ...]
n_Mc = n + n_w;
c_Mc = [c; zeros(n_w, 1)];          %%%%%I think no need to do this
H_Mc = diag([diag(H); zeros(n_w, 1)]);    %%%%I think no need to do this
lb_Mc = [lb; -1e10 * ones(n_w, 1)]; % No boundaries for w
ub_Mc = [ub; 1e10 * ones(n_w, 1)];
% ECs
Q_Mc = zeros(m, n_Mc);
for k = 1:m
    Qi = Q(n*(k-1)+1: n*k, :);
    loc = find(triu(Qi, 1)); % Nonzero index of Qi
    for j = 1 : length(loc)
        j_w = find(uidx_w == loc(j)); % Index corresponds to the j_wth value in w array
        Q_Mc(k, n + j_w) = 2 * Qi(loc(j)); % Set value for w
    end
end
A_Mc = [A, zeros(m, n_w)];
Aec_Mc = Q_Mc + A_Mc; % Matrix for equality constraints
bec_Mc = b;           % Vector for equality constraints
% ICs
Bic_Mc = zeros(4*n_w, n_Mc); % Matrix for inequality constraints
bic_Mc = zeros(4*n_w, 1);    % Vector for inequality constraints
for k = 1 : n_w
    loc = uidx_w(k);
    sz = [n n];
    [i, j] = ind2sub(sz, loc);
    xiL = lb(i); xjL = lb(j); xiU = ub(i); xjU = ub(j);
    Bic_Mc(4*(k-1)+1, [j,i,n+k]) = [xiL, xjL, -1]; bic_Mc(4*(k-1)+1) = xiL * xjL;
    Bic_Mc(4*(k-1)+2, [j,i,n+k]) = [xiU, xjU, -1]; bic_Mc(4*(k-1)+2) = xiU * xjU;
    Bic_Mc(4*(k-1)+3, [j,i,n+k]) = [-xiU,-xjL, 1]; bic_Mc(4*(k-1)+3) = -xiU * xjL;
    Bic_Mc(4*(k-1)+4, [j,i,n+k]) = [-xiL,-xjU, 1]; bic_Mc(4*(k-1)+4) = -xiL * xjU;
end

% 4. Compose the convex QP and solve it by quadprog function.
% x = quadprog(H,f,A,b,Aeq,beq,lb,ub)
x_lb = quadprog(2 * H_Mc,c_Mc, Bic_Mc, bic_Mc, Aec_Mc, bec_Mc, lb_Mc, ub_Mc);
f_lb = x_lb' * H_Mc * x_lb + c_Mc' * x_lb;

% 5. Apply fmincon to the original problem for an upper bound
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
fun = @(x) (x' * H * x + c' * x);
Bic_fm = []; bic_fm = []; Aec_fm = []; bec_fm = [];
x0 = zeros(n, 1);
nonlcon = @(x) original_ECs(n, m, Q, A, b, x);
x_ub = fmincon(fun, x0, Bic_fm, bic_fm, Aec_fm, bec_fm, lb, ub, nonlcon);
f_ub = fun(x_ub);

end