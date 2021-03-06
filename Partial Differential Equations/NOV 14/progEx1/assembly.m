function [A,B] = assembly(x,h,f)
    N = size(x,2)-2;
    A = zeros(N);
    for i= 1:N
        for j = 1:N
            if i==j
                A(i,j)= 1/h(i+1)+1/h(i);
            elseif i==j+1
                A(i,j)=-1/h(i);
            elseif i==j-1
                A(i,j)=-1/h(i+1);
            else
                A(i,j)=0;
            end
        end
    end
    B=zeros(N,1);
    if f ==1
        for i = 1:N
            B(i)=(x(i+2)-x(i))/2;
        end
    else
        for i = 1:N
            if -x(i+1)+0.5<h(i+1) && -x(i+1)+0.5>0
                a = 2*(x(i+2)-0.5)/h(i+1);
            else
                a=0;
            end
            if x(i+1)-0.5<h(i) && x(i+1)-0.5>0
                b = 2*(0.5-x(i))/h(i);
            else
                b=0;
            end
            B(i) = a+b;
            if x(i+1)==0.5
                B(i)= 1;
            end
        end


    end
end
