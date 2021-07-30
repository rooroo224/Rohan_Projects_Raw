#include "postProcessor.h"

/***************************************************************************************************
preProcessorControl
**************************************************************************************************/

void postProcessor::postProcessorControl(inputSettings* argSettings, triMesh* argMesh)
{
    mesh = argMesh;
    settings = argSettings;

    vtkVisualization();

    return;
}

/***************************************************************************************************
// Main visualization function
***************************************************************************************************/
void postProcessor::vtkVisualization()
{
    int nn = mesh->getNn();
    int ne = mesh->getNe();


     string dummy;

     // VTK Double Array
     vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
     pcoords->SetNumberOfComponents(nen);
     pcoords->SetNumberOfTuples(nn);

     //vtkDoubleArray type pcoords is filled with the data in meshPoints.
     for (int i=0; i<nn; i++)
         pcoords->SetTuple3(i,mesh->getNode(i)->getX(),mesh->getNode(i)->getY(),0.0f);

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

     // Previously collected data which are outputPoints, outputCells, scalarProperty, are written to
     // vtkUnstructuredGrid type grid.
     vtkSmartPointer<vtkUnstructuredGrid> unsGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
     unsGrid->SetPoints(outputPoints);
     unsGrid->SetCells(5,connectivity);

     //Whatever collected in unstructured grid above is written to the "Title_mype.vtu" file below.
     vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
     dummy = settings->getTitle();
     dummy.append(".vtu");
     cout << dummy << endl;
     writer->SetFileName(dummy.c_str());
     writer->SetInputData(unsGrid);
     writer->Write();

     // Now we write the "Title.pvtu" file which contains the informtaiton about other files.
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
