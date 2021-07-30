function [max] = inifnormcalc(u)

N = size(u,1);
max = 0;

for i=1:N
    for j=1:N
        if abs(u(i,j)) > max
            max = abs(u(i,j));
        end
    end
end

end
