meshOps = MeshOperations('unitSquare2ndOrder.msh');
nodes = meshOps.getNumberNodes();
lines = meshOps.getNumberOfLines();
tri = meshOps.getNumberOfTriangles();

[weight,points,ip]=meshOps.IntegrationRuleOfTriangle();
A = zeros(6,6,tri);
B = zeros(6,tri);
f = zeros(ip,1);

%%%creating local matrices
for i=1:1:tri
    detJ = meshOps.calcJacobianDeterminantOfTriangle(i);
    invJ = meshOps.calcInverseJacobianOfTriangle(i);
    shape = zeros(6,ip);
    delshape=zeros(ip,2,6);
    for j = 1:1:ip
        coord = meshOps.calcMappedIntegrationPointOfTriangle(i,[0,0])+meshOps.calcJacobianOfTriangle(i)*(points(j,:).');
        shape(1,j) = (1 - points(j,1) - points(j,2))*(1-2*points(j,1)-2*points(j,2));
        shape(2,j) = points(j,1)*(2*points(j,1)-1);
        shape(3,j) = points(j,2)*(2*points(j,2)-1);
        shape(4,j) = 4*points(j,1)*(1-points(j,1)-points(j,2));
        shape(5,j) = 4*points(j,1)*points(j,2);
        shape(6,j) = 4*points(j,2)*(1-points(j,1)-points(j,2));
        delshape(j,:,1) = [-3+4*points(j,1)+4*points(j,2) -3+4*points(j,1)+4*points(j,2)];
        delshape(j,:,2) = [4*points(j,1)-1 0];
        delshape(j,:,3) = [0 4*points(j,2)-1];
        delshape(j,:,4) = [4-8*points(j,1)-4*points(j,2) -4*points(j,1)];
        delshape(j,:,5) = [4*points(j,2) 4*points(j,1)];
        delshape(j,:,6) = [-4*points(j,2) 4-4*points(j,1)-8*points(j,2)];
        f(j) = sin(2*pi*coord(1));
    end
    for a=1:1:6
        for b=1:1:6
            for n=1:1:ip
                A(a,b,i)=A(a,b,i)+(weight(n)*dot((invJ*(delshape(n,:,a).')),(invJ*(delshape(n,:,b).')))*detJ);
            end
        end
    end
    for b=1:1:6
        for n=1:1:ip
            B(b,i)=B(b,i)+weight(n)*f(n)*shape(b,n)*detJ;
        end
    end
end
%%%%%

%%%assembling global matrices
K = zeros(meshOps.getNumberNodes()); %%global stiffness matrix
R = zeros(meshOps.getNumberNodes(),1); %global B matrix
for i=1:1:tri
    localset = meshOps.getNodeNumbersOfTriangle(i,2);
    K(localset,localset)= K(localset,localset)+ A(:,:,i);
    R(localset)= R(localset)+ B(:,i);
end
%%%

%%dirichlet boundary node list
nodevec = meshOps.getNodeNumbersOfLine(1,2);
for i=2:1:lines
    tag = meshOps.getTagOfLine(i);
    if(tag==2)
        nodevec = [nodevec,meshOps.getNodeNumbersOfLine(i,2)];
    end
end
%%%%

%%updating boundary values in global matrices
for j=1:1:nodes
    for i=1:1:length(nodevec)
            if (j==nodevec(i))
                R(j,1)=0;
                for b=1:1:nodes
                        if (b==j)
                            K(j,b)=1;
                        else K(j,b)=0;
                        end
                end
            end
    end
end
%%%%%%

u = meshOps.solve(K,R);
L2err = meshOps.L2errorOfPoissonProblem(u,2);
meshOps.plot(u,1)