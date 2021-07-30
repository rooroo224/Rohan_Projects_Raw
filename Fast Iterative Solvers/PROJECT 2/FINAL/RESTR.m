function [fnew] = RESTR(r, Nc)
% Fine to coarse
fnew=zeros(Nc+1,Nc+1);
for i=2:Nc
    ii=2*i-1;
    for j=2:Nc
        jj=2*j-1;
        fnew(i,j)=1/16*(r(ii-1,jj-1)+2*r(ii,jj-1)+r(ii+1,jj-1)+2*r(ii-1,jj)+4*r(ii,jj)+2*r(ii+1,jj)+r(ii-1,jj+1)+2*r(ii,jj+1)+r(ii+1,jj+1));
    end
end

end

