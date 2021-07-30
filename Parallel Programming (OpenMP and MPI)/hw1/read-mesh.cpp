/***************************************************************************************************
Name        : 2D_Unsteady_Diffusion.cpp
Author      : A. Emre Ongut
Version     : 2.1 - parallel
Copyright   : See the copyright notice in the README file.
Description : This code is designed for educational purposes.
             - For now, everything is for 2D linear triangular elements
             - uses Gauss quadrature integration with 7 points
             - supports constant Drichlet, Neumann and Robin boundary conditions
             - boundary conditions for ONLY FOR THE RECTANGULAR DOMAIN!
             - mixed boundary condition is not implemented!
***************************************************************************************************/

#include "settings.h"
#include "tri.h"
#include "postProcessor.h"

int main(int argc, char **argv) {
/***************************************************************************************************
MAIN PROGRAM FLOW
1. Pre-Processing Stage
    1.1. Settings
    1.2. Mesh
3. Post-Processing Stage
***************************************************************************************************/

    inputSettings*  settings    = new inputSettings(argc, argv);
    triMesh*        mesh        = new triMesh;
    postProcessor*  postP       = new postProcessor;

    // Pre-Processing Stage
    settings->prepareSettings();
    mesh->prepareMesh(settings);
//
    // Post-Processing Stage
    postP->postProcessorControl(settings, mesh);

    // Cleanup
    delete settings;
    delete mesh;
    delete postP;

    return 0;
}