function [etildeh] = PROLONG(etildeH,Nc)
% Coarse to fine

N = 2 * Nc;
etildeh=zeros(N + 1, N + 1);
for i = 2 : Nc
    ii = 2 * i - 1;
    for j = 2 : Nc
        jj = 2 * j - 1;
        etildeh(ii-1,jj-1) = etildeh(ii-1,jj-1) + (1/4*etildeH(i,j));
        etildeh(ii,jj-1)   = etildeh(ii,jj-1)   + (1/2*etildeH(i,j));
        etildeh(ii+1,jj-1) = etildeh(ii+1,jj-1) + (1/4*etildeH(i,j));
        etildeh(ii-1,jj)   = etildeh(ii-1,jj)   + (1/2*etildeH(i,j));
        etildeh(ii,jj)     = etildeh(ii,jj)     + (1*etildeH(i,j));
        etildeh(ii+1,jj)   = etildeh(ii+1,jj)   + (1/2*etildeH(i,j));
        etildeh(ii-1,jj+1) = etildeh(ii-1,jj+1) + (1/4*etildeH(i,j));
        etildeh(ii,jj+1)   = etildeh(ii,jj+1)   + (1/2*etildeH(i,j));
        etildeh(ii+1,jj+1) = etildeh(ii+1,jj+1) + (1/4*etildeH(i,j));     
    end
end

end

