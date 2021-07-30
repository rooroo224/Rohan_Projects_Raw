% Sijie Luo            - 403598 :)
% Rohan Krishna Balaji - 403596

function [val] = ano_case1()
% P1, P2, PP, EP,Power, Fuel, C, I1, I2, HE1, HE2, LE1, LE2, BF1, BF2, HPS, MPS, LPS
%  1   2   3   4    5     6   7   8   9   10   11   12   13   14   15   16   17   18 (INDEX)
%x = linprog(J,A,B,Aeq,Beq,lb,ub);
   
% Case 1: EP <= 12000
% Inequality Constraints
A = zeros(27,18); 
B = zeros(1,27);
for i=1:18
    A(i,i) = -1;            % making sure all variables are positive
end

A(19,8)  = 1;  A(20,7) = 1;   A(21,8) = 1;  A(21,10) = -1; A(22,9) = 1;  A(23,13) = 1;
A(24,17) = 1;  A(25,18) = -1; A(26,1) = -1; A(26,2) = -1;  A(26,4) = -1; A(27,4) = 1; 

B(19) = 87000;  B(20) = 28000;  B(21) = 60000; B(22) = 110000; B(23) = 64000;  B(24) = 123000;
B(25) = -45000; B(26) = -24500; B(27) = 12000;

% P1, P2, PP, EP,Power, Fuel, C, I1, I2, HE1, HE2, LE1, LE2, BF1, BF2, HPS, MPS, LPS
%  1   2   3   4    5     6   7   8   9   10   11   12   13   14   15   16   17   18 (INDEX)

%Mass Balance as equality constraints

Aeq = zeros(10,18); % We shall reshape the matrix into row matrix later
Beq = zeros(1,10);
Aeq(1,16) = 1; Aeq(1,8) = -1;  Aeq(1,9) = -1;  Aeq(1,14) = -1;
Aeq(2,8)  = 1; Aeq(2,10) = -1; Aeq(2,12) = -1; Aeq(2,7) = -1;
Aeq(3,9)  = 1; Aeq(3,11) = -1; Aeq(3,13) = -1;
Aeq(4,10) = 1; Aeq(4,11) = 1;  Aeq(4,14) = 1;  Aeq(4,15) = -1; Aeq(4,17) = -1;
Aeq(5,18) = 1; Aeq(5,12) = -1; Aeq(5,13) = -1; Aeq(5,15) = -1;

%Energy Balance as Equality Constriants
% P1, P2, PP, EP,Power, Fuel, C, I1, I2, HE1, HE2, LE1, LE2, BF1, BF2, HPS, MPS, LPS
%  1   2   3   4    5     6   7   8   9   10   11   12   13   14   15   16   17   18 (INDEX)

Aeq(6,8)   = 3163;   Aeq(6,10) = -2949;  Aeq(6,12) = -2911; Aeq(6,7) = -449; Aeq(6,1) = -3600;
Aeq(7,9)   = 3163;   Aeq(7,11) = -2949;  Aeq(7,13) = -2911; Aeq(7,2) = -3600;
Aeq(8,3)   = 1;      Aeq(8,1)  = -1;     Aeq(8,2)  = -1;
Aeq(9,5)   = 1;      Aeq(9,3)  = -1;     Aeq(9,4)  = -1;
Aeq(10,16) = 3163;   Aeq(10,6) = -0.75;

% Bounded Constraints
% P1, P2, PP, EP,Power, Fuel, C, I1, I2, HE1, HE2, LE1, LE2, BF1, BF2, HPS, MPS, LPS
%  1   2   3   4    5     6   7   8   9   10   11   12   13   14   15   16   17   18 (INDEX)

lb = []; ub = [];
lb(1,1) = 2500; lb(1,2) = 3000;
ub(1,1) = 6250; ub(1,2) = 9000;

% Optimization Objective Function
J = zeros(1,18);
J(1,6) = 1.5 * 10^-6; J(1,7) = 0.008; J(1,3) = 0.02; J(1,4) = 0.05-0.001;  

x = linprog(((0.001*12000)+J),A,B,Aeq,Beq,lb,ub);
val = (0.001*12000) + J*x;

end