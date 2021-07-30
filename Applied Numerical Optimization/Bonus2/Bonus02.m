%% Problem 1
clc; clear all;

n = 3;
m = 2;
c = [1; 1; 0];
b = [8;15];
A = zeros(2,3);
A(1,3) = 1;
H = zeros(3,3);
H(3,3) = 1;
Q1 = zeros(3,3);
Q1(1,2) = 0.5;
Q1(2,1) = 0.5;
Q2 = zeros(3,3);
Q2(3,2) = 0.5;
Q2(2,3) = 0.5;
Q = [Q1; Q2];
% bounds
lb = [0;0;0];
ub = [10;10;10];
[f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)

%% Problem 2
clc; clear all;

n = 4;
m = 3;
c = [1; 1; 0; 0];
b = [2; 3; 5];
A = zeros(3,4);
A(2,4) = 1;
A(3,1) = 1;
H = zeros(4,4);
H(3,3) = 1;
H(4,4) = 1;
Q1 = zeros(4,4);
Q1(1,2) = 0.5;
Q1(2,1) = 0.5;
Q1(2,3) = 0.5;
Q1(3,2) = 0.5;
Q2 = zeros(4,4);
Q2(1,2) = 0.5;
Q2(2,1) = 0.5;
Q3 = zeros(4,4);
Q3(2,3) = 0.5;
Q3(3,2) = 0.5;
Q = [Q1; Q2; Q3];
% bounds
lb = [0; 0; 0; 0];
ub = [10; 4; 10; 10];
[f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)

