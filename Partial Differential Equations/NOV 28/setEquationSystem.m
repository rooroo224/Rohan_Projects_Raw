function [K,R]=setEquationSystem(meshOps,fun)

tri = meshOps.getNumberOfTriangles();

[weight,points,ip]=meshOps.IntegrationRuleOfTriangle();
A = zeros(3,3,tri);
B = zeros(3,tri);
f = zeros(ip,1);
delshape=zeros(2,1,3);
delshape(:,:,1) = [-1;-1];
delshape(:,:,2) = [1;0];
delshape(:,:,3) = [0;1];

%%creating local element matrices
for i=1:1:tri
    detJ = meshOps.calcJacobianDeterminantOfTriangle(i);
    invJ = meshOps.calcInverseJacobianOfTriangle(i);
    shape = zeros(3,ip);
    
    for j = 1:1:ip
        coord = meshOps.calcMappedIntegrationPointOfTriangle(i,[0,0])+meshOps.calcJacobianOfTriangle(i)*(points(j,:).');
        shape(1,j) = 1 - points(j,1) - points(j,2);
        shape(2,j) = points(j,1);
        shape(3,j) = points(j,2);
        f(j) = fun(coord(1),coord(2));
    end
    for a=1:1:3
        for b=1:1:3
            for n=1:1:ip
                A(a,b,i)=A(a,b,i)+(weight(n)*dot((invJ*delshape(:,:,a)),(invJ*delshape(:,:,b)))*detJ);
            end
        end
    end
    for b=1:1:3
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
    localset = meshOps.getNodeNumbersOfTriangle(i,1);
    K(localset,localset)= K(localset,localset)+ A(:,:,i);
    R(localset)= R(localset)+ B(:,i);
end
