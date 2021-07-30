clear all; 
clc;
tic
%Basic Input Data 
n=4;l=4;nu1=2; nu2=1; gamma=2;
N=(2^n); h=1/N;

uh=zeros(N,N);ua=zeros(N,N);
%defining the RHS function
 f=zeros(N,N);
 for j=2:N
     for i=2:N
         f(i,j)=8*pi^2*sin(2*pi*h*i)*sin(2*pi*h*j);
         ua(i,j)=sin(2*pi*i*h)*sin(2*pi*j*h);
     end
 end
















