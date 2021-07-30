function [A,B]=setBoundaryValues(meshOps,K,R,nodevec,g_d)

nodes = meshOps.getNumberNodes();
A = K;
B = R;

for j=1:1:nodes
    for i=1:1:length(nodevec)
            if (j==nodevec(i))
                B(j,1)=g_d;
                for b=1:1:nodes
                        if (b==j)
                            A(j,b)=1;
                        else A(j,b)=0;
                        end
                end
            end
    end
end
