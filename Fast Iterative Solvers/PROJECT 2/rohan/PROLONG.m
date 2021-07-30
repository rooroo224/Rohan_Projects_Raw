%function [eh] = PROLONG(eH,Nc,eh)
function [eh] = PROLONG(eH,Nc,N)
    %NO FOR TEST
    eh=zeros(N,N);
    for i=2:Nc-1
        ii=2*i;
        for j=2:Nc-1
            jj=2*j;
            eh(ii-1,jj-1)= eh(ii-1,jj-1)+(1/4*eH(i,j));
            eh(ii,jj-1)  = eh(ii,jj-1)  +(1/2*eH(i,j));
            eh(ii+1,jj-1)= eh(ii+1,jj-1)+(1/4*eH(i,j));
            eh(ii-1,jj)  = eh(ii-1,jj)  +(1/2*eH(i,j));
            eh(ii,jj)    = eh(ii,jj)    +(1*eH(i,j));
            eh(ii+1,jj)  = eh(ii+1,jj)  +(1/2*eH(i,j));
            eh(ii-1,jj+1)= eh(ii-1,jj+1)+(1/4*eH(i,j));
            eh(ii,jj+1)  = eh(ii,jj+1)  +(1/2*eH(i,j));
            eh(ii+1,jj+1)= eh(ii+1,jj+1)+(1/4*eH(i,j));     
        end
    end
end
