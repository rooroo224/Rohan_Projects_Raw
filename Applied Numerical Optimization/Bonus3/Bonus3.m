clc; close all; clear all;

% Sijie Luo - 403598
% Rohan Krishna Balaji - 403596

% USER DEFINED M (numbe of intervals)
% Task 1
m(1) = 10;
m(2) = 50;
m(3) = 100;

% Task 2
%Printing the objective functions
for i = 1:size(m,2)
    disp('Calculating...')
    fprintf('The optimal objective function values with M = %i is %f\n',m(i),definitions(m(i)))
end


% Task 3:

% In a single shooting (late discretization) approach, the objective value in general
% decreases with a finer discretization of the control variables.

% Reasoning:
% In single shoting method we use control variables as linear
% combination of basis functions, this appproximation gets more accurate
% with finer discretization. Since our objective function is minimization,
% on better approximation (beacause of fine discretization) the objecive
% function decreases.