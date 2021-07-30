function u = linEqnsolver(A,B)

    N = size(B,1);
    u(1)=0;
    u(N+2)=0;
    u(2:N+1) = A\B;
    
end