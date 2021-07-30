function [eh] = PROLONG(eH,Nc)
N=2*Nc; eh=zeros(N+1,N+1);
for i=2:Nc
    ii=2*i-1;
    for j=2:Nc
        jj=2*j-1;
        eh(ii-1,jj-1) = eh(ii-1,jj-1) + (1/4*eH(i,j));
        eh(ii,jj-1)   = eh(ii,jj-1)   + (1/2*eH(i,j));
        eh(ii+1,jj-1) = eh(ii+1,jj-1) + (1/4*eH(i,j));
        eh(ii-1,jj)   = eh(ii-1,jj)   + (1/2*eH(i,j));
        eh(ii,jj)     = eh(ii,jj)     + (1*eH(i,j));
        eh(ii+1,jj)   = eh(ii+1,jj)   + (1/2*eH(i,j));
        eh(ii-1,jj+1) = eh(ii-1,jj+1) + (1/4*eH(i,j));
        eh(ii,jj+1)   = eh(ii,jj+1)   + (1/2*eH(i,j));
        eh(ii+1,jj+1) = eh(ii+1,jj+1) + (1/4*eH(i,j));     
    end
end

end

