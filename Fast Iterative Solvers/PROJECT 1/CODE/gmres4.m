function [x,rho,res]=gmres4(values,colid,rowptr,rows,b,xo,m,count)
    [z,q]=size(b); ro=zeros(z,1); v=zeros(z,m); g=zeros(m+1,1); e1=zeros(m+1,1); e1(1,1)=1;
    c=zeros(m,1); s=zeros(m,1);temp=zeros(z,1);r=zeros(z,1); 
%calculating initial error ro based on initial vector xo
temp= mat_vect_prod2(values,colid,rowptr,rows,xo);
for k=1:z
    ro(k,1)=b(k,1)-temp(k,1);
end
normro=0;
for k=1:z
    normro=normro+ ro(k,1)*ro(k,1);
end
normro= sqrt(normro);
for k=1:z
    v(k,1)=ro(k,1)/normro;
end
%done with calculating v1
    g= normro*e1;h=zeros(m+1,m);
for j=1:m
      %GetKrylov part
    pass=zeros(z,1);
    for p=1:z
        pass(p,1)=v(p,j); %just coping vj into a vector so that mat-vect prod could be done
    end
    w= mat_vect_prod2(values,colid,rowptr,rows,pass);
    for i=1:j
        temp=0; 
         for p=1:z
            temp=temp+v(p,i)*w(p,1); %dot product to calculate hij
         end
         h(i,j)=temp; 
         for k=1:z
             w(k,1)=w(k,1)-(h(i,j)*v(k,i)); %updating w using new hij
         end
    end
       normw=0;
    for p=1:z
        normw=normw+w(p,1)*w(p,1); %calculating norm of w
    end
    h(j+1,j)=sqrt(normw);
    for k=1:z
        v(k,j+1)=w(k,1)/h(j+1,j);  %calculating new orthogonal vector vj+1
    end
    for k=2:j
        temp1= h(k-1,j);  
        h(k-1,j)= ((c(k-1,1)*h(k-1,j)) + (s(k-1,1)*h(k,j))) ;
        h(k,j)= ((-s(k-1,1)*temp1)+ (c(k-1,1)*h(k,j)));
    end
    c(j,1)= h(j,j)/sqrt(h(j,j)^2+h(j+1,j)^2);
    s(j,1)= h(j+1,j)/sqrt(h(j,j)^2+h(j+1,j)^2);
    h(j,j)=c(j,1)*h(j,j)+s(j,1)*h(j+1,j);
    g(j+1,1)= -s(j,1)*g(j,1);
    g(j,1)=c(j,1)*g(j,1);
    res(j+(count*m),1)=abs(g(j+1,1))/normro;
end
 y=zeros(m,1); 
    y(m,1)=g(m)/h(m,m); 
   for k=m-1:-1:1
       dummy=0; 
          for j=k+1:m
              dummy =dummy+h(k,j)*y(j,1);
          end
       y(k,1)=(g(k,1)-dummy)/h(k,k);
   end
   V=zeros(z,m);
   for i=1:z
       for j=1:m
           V(i,j)=v(i,j);
       end
   end
  x= xo+ mat_vect_prod(V,y);
 % finding rho to restart gmres
  temp= mat_vect_prod2(values,colid,rowptr,rows,x);
for k=1:z
    r(k,1)=b(k,1)-temp(k,1);
end
rho=0;
for k=1:z
     rho= rho+ r(k,1)*r(k,1);
end
rho= sqrt(rho);











