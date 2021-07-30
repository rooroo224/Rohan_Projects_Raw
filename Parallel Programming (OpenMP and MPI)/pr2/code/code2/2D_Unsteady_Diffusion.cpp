/***************************************************************************************************
 * Name        : 2D_Unsteady_Diffusion.cpp
 * Description : This code is designed for educational purposes.
 *            - For now, everything is for 2D linear triangular elements
 *            - uses Gauss quadrature integration with 7 points
 *            - supports constant Drichlet
 ***************************************************************************************************/

#include "settings.h"
#include "tri.h"
#include "solver.h"
#include "postProcessor.h"
#include "omp.h"
#include <ctime>

int main(int argc, char **argv) {
/******************************************************************************* ********************
 * MAIN PROGRAM FLOW
 * 1. Pre-Processing Stage
 *    1.1. Settings
 *    1.2. Mesh
 * 2. Solution Stage
 * 3. Post-Processing Stage
 ***************************************************************************************************/

    double starttime = omp_get_wtime();

    inputSettings*  settings    = new inputSettings(argc, argv);
    triMesh*        mesh        = new triMesh;
    femSolver*      solver      = new femSolver;
    postProcessor*  postP       = new postProcessor;

    // Pre-Processing Stage
    settings->prepareSettings();
    mesh->prepareMesh(settings);

    // Solution Stage
    solver->solverControl(settings, mesh);

    // Post-Processing Stage
    postP->postProcessorControl(settings, mesh);

    double endtime = omp_get_wtime();
    cout << "Elapsed time is " << fixed << double(endtime - starttime) << endl;

    // Cleanup
    delete settings;
    delete mesh;
    delete solver;
    delete postP;

    return 0;
}
