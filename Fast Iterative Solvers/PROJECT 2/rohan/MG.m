function [uout,count]= MG(l,ul, fl, gamma, nu1, nu2, N,count)
    count=count+1;
    Nc=N/2;
    e=zeros(Nc,Nc); h=1/N;
    ul= GS(ul,fl,nu1,N); 
    r=zeros(N,N);
    
    %Calculating the Residual Before Restricion operator
    for i=2:N-1
        for j=2:N-1
            r(i,j)=fl(i,j)+ ((ul(i-1,j)-2*ul(i,j)+ul(i+1,j))/h^2) + ((ul(i,j-1)-2*ul(i,j)+ul(i,j+1))/h^2);
        end
    end
    global a
    a(count,1)=norm(r,inf);
    [rl] = RESTR(r,Nc);
    if l==1
        [e]= GS(e,rl,nu1,Nc);
    else
        e=zeros(Nc,Nc);
        for j=1:gamma
            [e,count]= MG(l-1,e,-rl,gamma, nu1,nu2,Nc,count);
        end
    end
    e= PROLONG(e,Nc,N);
    uout = GS(ul-e,fl,nu2,N);
end








