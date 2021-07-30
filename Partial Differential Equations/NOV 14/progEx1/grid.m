function [x,h]=grid(N,L)
 
 %  x = linspace(0,L,N+2); %For equal spacing
 
 %%%%%%%%%%%%%%%
 % for unequal spacing
 
    x(1)=0;
    x(N+2)=L;
    x(2:N+1) = sort(rand(1,N)); 
    
 %%%%%%%%%%%%%%%
 
 h = zeros(1,N+1);
    for i=1:N+1
        h(i)=x(i+1)-x(i); %grid size
    end
end