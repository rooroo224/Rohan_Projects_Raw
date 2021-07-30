function [rH] = RESTR(rh, Nc)
% Fine to coarse
rH=zeros(Nc+1,Nc+1);
rH = zeros(Nc + 1);
for i = 2 : Nc
    ii = 2 * i - 1;
    for j = 2 : Nc
        jj = 2 * j - 1;
        rH(i,j) = 1/16 * (rh(ii-1,jj-1) + 2 * rh(ii,jj-1) + rh(ii+1,jj-1) + ...
                        2 * rh(ii-1,jj) + 4 * rh(ii,jj) + 2 * rh(ii+1,jj) + ...
                        rh(ii-1,jj+1) + 2 * rh(ii,jj+1) + rh(ii+1,jj+1));
    end
end

end

