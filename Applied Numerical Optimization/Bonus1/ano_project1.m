% Sijie Luo            - 403598 :)
% Rohan Krishna Balaji - 403596

clear; close all; clc;

val1 = ano_case1();
val2 = ano_case2();

if (val1>val2)
    fprintf('Case 2 is better = %f\n', val2);
else
    fprintf('Case 1 is better = %f\n ', val1);
end
