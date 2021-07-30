#include "solver.h"

/***************************************************************************************************
solverControl
****************************************************************************************************
Flow control to prepare and solve the system.
***************************************************************************************************/
void femSolver::solverControl(inputSettings* argSettings, triMesh* argMesh)
{
    cout << endl << "=================== SOLUTION =====================" << endl;

    mesh = argMesh;
    settings = argSettings;
    int ne = mesh->getNe();

    for(int iec=0; iec<ne; iec++)
    {
        calculateJacobian(iec);
        calculateElementMatrices(iec);
    }
    applyDrichletBC();
    accumulateMass();
    explicitSolver();

    return;
}

/***************************************************************************************************
calculateJacobian
****************************************************************************************************
Compute and store the jacobian for each element.
***************************************************************************************************/
void femSolver::calculateJacobian(const int e)
{
    int myNode;     // node number for the current node
    double x[nen];  // x values for all the nodes of an element
    double y[nen];  // y values for all the nodes of an element

    triMasterElement* ME = mesh->getME(0);  // for easy access to the master element
    double * xyz = mesh->getXyz();

    double J[2][2];     // Jacobian for the current element
    double detJ;        // Jacobian determinant for the current element
    double invJ[2][2];  // inverse of Jacobian for the current element
    double dSdX[3];     // dSdx on a GQ point
    double dSdY[3];     // dSdy on a GQ point

    // collect element node coordinates in x[3] and y[3] matrices
    for (int i=0; i<nen; i++)
    {
        myNode =  mesh->getElem(e)->getConn(i);
        x[i] = xyz[myNode*nsd+xsd];
        y[i] = xyz[myNode*nsd+ysd];
    }

    // for all GQ points detJ, dSDx[3] and dSdY[3] are determined.
    for (int p=0; p<nGQP; p++)
    {
        // Calculate Jacobian
        J[0][0] = ME[p].getDSdKsi(0)*x[0] + ME[p].getDSdKsi(1)*x[1] + ME[p].getDSdKsi(2)*x[2];
        J[0][1] = ME[p].getDSdKsi(0)*y[0] + ME[p].getDSdKsi(1)*y[1] + ME[p].getDSdKsi(2)*y[2];
        J[1][0] = ME[p].getDSdEta(0)*x[0] + ME[p].getDSdEta(1)*x[1] + ME[p].getDSdEta(2)*x[2];
        J[1][1] = ME[p].getDSdEta(0)*y[0] + ME[p].getDSdEta(1)*y[1] + ME[p].getDSdEta(2)*y[2];

        //Calculate determinant of Jacobian and store in mesh
        detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
        mesh->getElem(e)->setDetJ(p, detJ);

        // Calculate inverse of Jacobian
        invJ[0][0] =  J[1][1]/detJ;
        invJ[0][1] = -J[0][1]/detJ;
        invJ[1][0] = -J[1][0]/detJ;
        invJ[1][1] =  J[0][0]/detJ;

        // Calculate dSdx and dSdy and store in mesh
        dSdX[0] = invJ[0][0]*ME[p].getDSdKsi(0) + invJ[0][1]*ME[p].getDSdEta(0);
        dSdX[1] = invJ[0][0]*ME[p].getDSdKsi(1) + invJ[0][1]*ME[p].getDSdEta(1);
        dSdX[2] = invJ[0][0]*ME[p].getDSdKsi(2) + invJ[0][1]*ME[p].getDSdEta(2);
        dSdY[0] = invJ[1][0]*ME[p].getDSdKsi(0) + invJ[1][1]*ME[p].getDSdEta(0);
        dSdY[1] = invJ[1][0]*ME[p].getDSdKsi(1) + invJ[1][1]*ME[p].getDSdEta(1);
        dSdY[2] = invJ[1][0]*ME[p].getDSdKsi(2) + invJ[1][1]*ME[p].getDSdEta(2);

        mesh->getElem(e)->setDSdX(p, 0, dSdX[0]);
        mesh->getElem(e)->setDSdX(p, 1, dSdX[1]);
        mesh->getElem(e)->setDSdX(p, 2, dSdX[2]);
        mesh->getElem(e)->setDSdY(p, 0, dSdY[0]);
        mesh->getElem(e)->setDSdY(p, 1, dSdY[1]);
        mesh->getElem(e)->setDSdY(p, 2, dSdY[2]);
    }

    return;
}

/***************************************************************************************************
void femSolver::calculateElementMatrices(const int e)
****************************************************************************************************
Compute the K, M and F matrices. Then accumulate the total mass into the node structure.
***************************************************************************************************/
void femSolver::calculateElementMatrices(const int e)
{
    // Add code here
  double M[nen][nen] = {0.0};
  double K[nen][nen] = {0.0};
  double F[nen] = {0.0};
  double a = settings->getD();     //thermal diffusivity
  double * xyz = mesh->getXyz();
  double f = settings->getSource();     // source  term

 
  // compute element level matrices for K,M,F
  for(int k=0; k < nGQP; k++)     // loop over quadrature points
    {
      double w = mesh->getME(k)->getWeight();
      double J = mesh->getElem(e)->getDetJ(k);

      for(int i=0; i<nen; i++)
	{
	  double Si = mesh->getME(k)->getS(i);
	  double Six = mesh->getElem(e)->getDSdX(k,i);
	  double Siy = mesh->getElem(e)->getDSdY(k,i);

	  for(int j=0; j<nen; j++)
	    {
	      double Sj = mesh->getME(k)->getS(j);
	      double Sjx = mesh->getElem(e)->getDSdX(k,j);
	      double Sjy = mesh->getElem(e)->getDSdY(k,j);
	      M[i][j] = M[i][j] + w*(Si*Sj)*J;
	      K[i][j] = K[i][j] + a*w*(Six*Sjx + Siy*Sjy)*J;
	    }

	  int n = mesh->getElem(e)->getConn(i);     //global node number
	  double x = xyz[n*nsd + xsd];
	  double y = xyz[n*nsd + ysd];
	  double r = sqrt(x*x + y*y);
	  if(r<0.01)     // check if node distance is less than or equal to R1 with 1% error
	    {
	      F[i] = F[i] + w*Si*J*f;
	    }
	  else
	    {
	      F[i]=0;
	    }
	}
    }

  double Me = 0;     //sum of M matrix
  double Mj = 0;     // diagonal sum of M matrix
  double ML[nen];     //lumped mass vector

  for(int i=0; i<nen; i++)
    {
      mesh->getElem(e)->setF(i,F[i]);     //set F vector
      Mj = Mj + M[i][i];
      for(int j=0; j<nen; j++)
	{
	  Me = Me + M[i][j];
	  mesh->getElem(e)->setK(i,j,K[i][j]);     //set K matrix
	}
    }

  //calculate lumped mass vector
  for(int i=0; i<nen; i++)
    {
      int n = mesh->getElem(e)->getConn(i);     //global node number
      ML[i] = M[i][i]*(Me/Mj);
      mesh->getNode(n)->addMass(ML[i]);     //add lumped mass for node n
      mesh->getElem(e)->setM(i,ML[i]);     //set M matrix

    }

  if(e==3921)
    {
    
      cout<<"K: "<<endl;
      for(int i=0; i<nen; i++)
	{
	  for(int j=0; j<nen; j++)
	    {
	      cout<<K[i][j]<<" ";
	    }
	  cout<<endl;
	}
      cout<<"M: "<<endl;
      for(int i=0;i<nen;i++)
	cout<<ML[i]<<endl;
      cout<<"F: "<<endl;
      for(int i=0; i<nen;i++)
	cout<<F[i]<<endl;
    }


    return;
}

/***************************************************************************************************
void femSolver::applyDrichletBC()
****************************************************************************************************
if any of the boundary conditions set to Drichlet type
    visits all partition level nodes
        determines if any the nodes is on any of the side surfaces

***************************************************************************************************/
void femSolver::applyDrichletBC()
{
    // Add code here

    int nn = mesh->getNn();
    double * xyz = mesh->getXyz();
    double * T = mesh->getT();     // global temperature array
    double x,y;
    int d=0;     // number of dirichlet nodes
    for(int i=0; i<nn; i++)
      {
	x = xyz[i*nsd+xsd];     // x and y coordinates for node 'i'
	y = xyz[i*nsd+ysd];
	double r = sqrt(x*x + y*y);
	double tol = (r - 0.1)*100/0.1;     //tolerance of 0.1%
	if(tol >= -0.1 && tol<= 0.1)
	  {
	    d = d+1;
	    mesh->getNode(i)->setBCtype(1);     // Dirichlet BC
	    T[i] = settings->getBC(1)->getValue1();
	  }
      }
    cout<<"Dirichlet Nodes = "<<d<<endl;
    return;
}

/***************************************************************************************************
* void femSolver::explicitSolver()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::explicitSolver()
{
    // Add code here

  int iter = settings->getNIter();
  double dt = settings->getDt();
  int ne = mesh->getNe();
  int nn = mesh->getNn();
  double * MTnew = mesh->getMTnew();
  double * T = mesh->getT();
  double * massG = mesh->getMassG();
  double scale,res,gerror;
  double lim = pow(10,-7);
  int i;

  for(i=1; i<=iter; i++)      //loop over iterations
    {

      for(int j=0; j<ne; j++)     // loop over Elements
	{
	  double * M = mesh->getElem(j)->getMptr();
	  double * K = mesh->getElem(j)->getKptr();
	  double * F = mesh->getElem(j)->getFptr();

	  for(int k=0; k<nen; k++)     // calculate MTnew
	    {
	      int n = mesh->getElem(j)->getConn(k);     //global node number
	      double KT = 0.0;
	      double MT = 0.0;
	      
	      for(int l=0; l<nen; l++)
		{
		  int p = mesh->getElem(j)->getConn(l);     //global node number
		  KT = KT + K[k*nen+l]*T[p];
		}

	      MT = M[k]*T[n];
	      MTnew[n] = MTnew[n]+ MT + dt*(F[k] - KT);   // MTnew for global node 'n'

	    }
	}

      gerror = 0.0;   // global error
      int nd = 0;     // non-Dirichlet nodes

      for(int j=0; j<nn; j++)     // loop over nodes
	{
	  int type = mesh->getNode(j)->getBCtype();
	  if(type!=1)     //non-dirichlet node
	    {
	      nd = nd +1;
	      double Tnew = MTnew[j]/massG[j];     // new Temp at current node
	      double error = (Tnew-T[j])*(Tnew-T[j]);
	      gerror = gerror + error;
	      T[j] = Tnew;     // store new Temperature
	    }

	  MTnew[j]=0.0;    //make 0 for next outer loop iteration
	}

      gerror = sqrt(gerror/nd);     // compute global error

      if(i==1)
	{
          scale = gerror;     // scaled value for 1st iteration
	  cout<<"Non-Dirichlet Nodes = "<<nd<<endl;
	}
      res = gerror/scale;     //rms error or residual
      if(res<lim) 
	break;
      //   cout<<"     "<<res<<"      "<<i<<endl;
    }

  cout<<"RMS Error = "<<res<<endl;
  cout<<"Iterations run = "<<i<<endl;

    return;
}


/***************************************************************************************************
* void femSolver::accumulateMass()
***************************************************************************************************
*
**************************************************************************************************/
void femSolver::accumulateMass()
{
    int nn = mesh->getNn();
    double * massG = mesh->getMassG();


    for(int i=0; i<nn; i++)
    {
        massG[i] = mesh->getNode(i)->getMass();
    }

    return;
}
