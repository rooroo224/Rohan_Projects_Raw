function [x,rho,count,e]= restart_gmres2(values,colid,rowptr,rows,xo,b,m)
[z,~]=size(b); ro=zeros(z,1); v=zeros(z,m); temp=zeros(z,1); x=zeros(z,1); xm=zeros(z,1);tol=10^-8;
xstar=ones(1030,1);
%calculating initial error ro based on initial vector xo
temp= mat_vect_prod2(values,colid,rowptr,rows,xo); count=0;
for z=1:z
    ro(z,1)=b(z,1)-temp(z,1);
end
rhoi=0;
for k=1:z
    rhoi=rhoi+ ro(k,1)*ro(k,1);
end
rhoi= sqrt(rhoi);
for i=1:z
    x(i,1)= xo(i,1);
rho=rhoi;    
end
while rho> tol
    [xm, rho]= gmres2(values,colid,rowptr,rows,b,x,m);
    rho=rho/rhoi
    count=count+1
    for k=1:z
        x(k,1)=xm(k,1);
        e(k,count)=x(k,1)-xstar(k,1);
    end
end
end



