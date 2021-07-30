function [uout,count] = MG(l, ul, fl, gamma, nu1, nu2, count)
count=count+1; [len,~]=size(ul); N=len-1; Nc=N/2; h=1/N;
r=zeros(size(ul)); e=zeros(Nc + 1);
ul=GS(ul, fl, nu1);

for i=2:N 
    for j=2:N
        r(i,j)=fl(i,j)+((ul(i-1,j)-2*ul(i,j)+ul(i+1,j))/(h^2))+((ul(i,j-1)-2*ul(i,j)+ul(i,j+1))/(h^2)) ;
    end
end
global a
a(count,1)=norm(r,inf);
rl=RESTR(r,Nc);  
if l==1
    e=GS(e,-rl,1);
else
    e=zeros(Nc + 1);  
    for i=1:gamma
        [e,count]=MG(l-1,e,-rl,gamma,nu1,nu2,count);
    end
end
e=PROLONG(e, Nc); 
uout=GS(ul-e,fl,nu2);
end

