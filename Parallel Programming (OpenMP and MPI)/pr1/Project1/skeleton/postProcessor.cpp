#include "postProcessor.h"
#include "math.h"

/***************************************************************************************************
preProcessorControl
**************************************************************************************************/

void postProcessor::postProcessorControl(inputSettings* argSettings, triMesh* argMesh)
{
    mesh = argMesh;
    settings = argSettings;

    cout << endl << "================ POST-PROCESSING =================" << endl;

    evaluateLimits();
    compareAnalytical();
    vtkVisualization();

    return;
}

/***************************************************************************************************
Evaluates the maximum and minimum temperatures in the field
***************************************************************************************************/
void postProcessor::evaluateLimits()
{
    int nn = mesh->getNn();
    double * T = mesh->getT();
    double Tmp;
    minT = std::numeric_limits<double>::max();
    maxT = std::numeric_limits<double>::min();

    for(int i=0; i<nn; i++)
    {
        Tmp = T[i];
        if(Tmp < minT)
            minT = Tmp;
        if(Tmp > maxT)
            maxT = Tmp;
    }
    cout << "Tmin " << minT << endl << "Tmax " << maxT << endl;

    return;
}

/***************************************************************************************************
Compare with analytical solution
***************************************************************************************************/
void postProcessor::compareAnalytical()
{
    int nn = mesh->getNn();
    double * T = mesh->getT();
    double * xyz = mesh->getXyz();
    double x, y, radius, Ta;
    // Add code here
    double To=300;double  Q=10; double alpha= settings->getD();
    double pi=3.14; double coeff;
    double error=0;  x=0; y=0; radius=0;     //declaring and initialising req variables
    coeff= (double)(Q/(2*pi*alpha));
    for(int i=0; i<nn;i++)
      {
	x= xyz[i*nsd+xsd];        //Loop over the nodes and finding analytical temperature
	y= xyz[i*nsd+ysd];        //based on distace from center
	radius= sqrt((x*x +y*y));
	Ta=0;
	if (radius< 0.01)
	  {
	    Ta= To-(coeff)*(0.5*(((double)(radius*radius)/0.01*0.01)-1)+log(0.1)); 
	  }
	else if (radius>=0.01)
	  {
	    Ta= To-(coeff)*log(radius/0.1);
	  }
	else
	  cout<<"error in code"<<endl;
	
	error=error+ ((T[i]-Ta)*(T[i]-Ta));  //Calculating disceretzation error wrt to 
      }                                      //Analytical and Numerical temperature
    error= sqrt((double)(error)/nn);
    cout<<" Discretization Error= "<<error<<endl;

    return;
}

/***************************************************************************************************
// Main visualization function
***************************************************************************************************/
void postProcessor::vtkVisualization()
{
    int nn = mesh->getNn();
    int ne = mesh->getNe();
    double * T = mesh->getT();
    double * xyz = mesh->getXyz();

    string dummy;

    // VTK Double Array
    vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
    pcoords->SetNumberOfComponents(nen);
    pcoords->SetNumberOfTuples(nn);

    //vtkDoubleArray type pcoords is filled with the data in meshPoints.
    for (int i=0; i<nn; i++)
        pcoords->SetTuple3(i,xyz[i*nsd+xsd],xyz[i*nsd+ysd],0.0f);

    //vtkPoints type outputPoints is filled with the data in pcoords.
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetData(pcoords);

    //Connectivity is written to vtkCellArray type outputCells
    vtkSmartPointer<vtkCellArray> connectivity = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i<ne; i++)
    {
        connectivity->InsertNextCell(nen);
        for(int j=0; j<nen; j++)
            connectivity->InsertCellPoint(mesh->getElem(i)->getConn(j));
    }

    // Scalar property
    vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetName("Temperature");
    for(int i=0; i<nn; i++)
        scalar->InsertNextValue(T[i]);

    // Previously collected data which are outputPoints, outputCells, scalarProperty, are written to
    // vtkUnstructuredGrid type grid.
    vtkSmartPointer<vtkUnstructuredGrid> unsGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unsGrid->SetPoints(outputPoints);
    unsGrid->SetCells(5,connectivity);
    unsGrid->GetPointData()->SetScalars(scalar);

    //Whatever collected in unstructured grid above is written to the "Title_mype.vtu" file below.
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    dummy = settings->getTitle();
    dummy.append(".vtu");
    cout << dummy << endl;
    writer->SetFileName(dummy.c_str());
    writer->SetInputData(unsGrid);
    writer->Write();

    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
    dummy = settings->getTitle();
    dummy.append(".pvtu");
    pwriter->SetFileName(dummy.c_str());
    #if VTK_MAJOR_VERSION <= 5
        pwriter->SetInput(unsGrid);
    #else
        pwriter->SetInputData(unsGrid);
    #endif
    pwriter->Write();

    return;
}
