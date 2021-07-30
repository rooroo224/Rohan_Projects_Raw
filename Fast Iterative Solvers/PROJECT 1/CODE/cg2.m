function [x,e,count]= cg2(values,colid,rowptr,rows,b,x)
    %function call and initialising required vectors
    [z,~]=size(b); r= zeros(z,1); p= zeros(z,1); temp=zeros(z,1); alpha=0; beta=0; vect=zeros(z,1);
     xstar=ones(5357,1);
    temp= mat_vect_prod3(values,colid,rowptr,rows,x);
    %calculating initial residual
    for k=1:z
        r(k,1)=b(k,1)-temp(k,1);
    end
     for k=1:z
        p(k,1)=r(k,1);
     end
    normro=0;
    %finding Norm
    for k=1:z
         normro=normro+ r(k,1)*r(k,1);
    end
    normro= sqrt(normro);
    normr=normro;
    count=0;
    %check for convergence 
    while normr>10^-8
        count=count+1
        dot1=0;
        for k=1:z
            dot1=dot1+r(k,1)*r(k,1);
        end
        vect= mat_vect_prod3(values,colid,rowptr,rows,p);
        dot2=0;
        for k=1:z
            dot2=dot2+vect(k,1)*p(k,1);
        end
        %finding alpha
        alpha= dot1/dot2;
       for k=1:z
           x(k,1)=x(k,1)+alpha*p(k,1);
           r(k,1)=r(k,1)-alpha*vect(k,1);
           e(k,count)=x(k,1)-xstar(k,1);
       end
       dot3=0;
        for k=1:z
            dot3=dot3+r(k,1)*r(k,1);
        end
        beta=dot3/dot1;
        for k=1:z
            p(k,1)=r(k,1)+beta*p(k,1);
        end
        normr=0;
        for k=1:z
            normr=normr+ r(k,1)*r(k,1);
        end
        normr= sqrt(normr);
        normr=normr/normro;
        res(count,1)=normr;
    end
    %plotting
    figure(1)
    [fg,~]=size(res)
    semilogy(1:fg,res(:))
    xlabel('Number of Iterations');
    ylabel('Standard 2 Norm (semilog)');
    title('Congugate Gradient-||r||_2');
end