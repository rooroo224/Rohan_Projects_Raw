meshOps = MeshOperations('unitSquare2.msh');
nodes = meshOps.getNumberNodes();
lines = meshOps.getNumberOfLines();
tri = meshOps.getNumberOfTriangles();
fun = @(x,y) sin(2*pi*x);

[K,R] = setEquationSystem(meshOps,fun); %create global stiffness and b matrices

%%list of dirichlet boundary nodes
nodevec = meshOps.getNodeNumbersOfLine(1,1);
for i=2:1:lines
    tag = meshOps.getTagOfLine(i);
    if(tag==2)
        nodevec = [nodevec,meshOps.getNodeNumbersOfLine(i,1)];
    end
end
%%%%

g_dir = 0; %%initialize dirichlet boundary value

[A,B] = setBoundaryValues(meshOps,K,R,nodevec,g_dir); %update boundary values in matrices
                
u = meshOps.solve(A,B); 
L2err = meshOps.L2errorOfPoissonProblem(u,1);
meshOps.plot(u,1)
    