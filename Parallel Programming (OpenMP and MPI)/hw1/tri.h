#ifndef TRI_H_
#define TRI_H_

#include "settings.h"

/***************************************************************************************************
NODE LEVEL DATA STRUCTURE
****************************************************************************************************
Each node has
    - coordinates
    - a variable
Comments:
    - For the current application, the only variable is temperature since we are solving transient
      diffusion equation.
***************************************************************************************************/
class triNode
{
    private:
        friend class femSolver;
        // VARIABLES
        double x; // x-coordinate
        double y; // y-coordinate
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        triNode(){x=0.0f; y=0.0f; };
        // SETTERS
        inline void   setX (double value) {x = value;};
        inline void   setY (double value) {y = value;};

        // GETTERS
        inline double getX() {return x;};
        inline double getY() {return y;};
};

/***************************************************************************************************
ELEMENT LEVEL DATA STRUCTURE
****************************************************************************************************
An element has
    - connectivity [nen]
    - boundary condition [nef]
    - Jacobian determinant [nGQP]
    - dsdX [3xnGQP]
    - dSdY [3xnGQP]
***************************************************************************************************/
class triElement
{
    private:
        // VARIABLES
        int conn[nen];
        int lConn[nen];
        int FG[nef];
        struct{
            double detJ;
            double dSdX[3];
            double dSdY[3];
        }GQP[nGQP];
        double M[nen];
        double F[nen];
        double K[nen][nen];
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        // SETTERS
        void setConn (int i, int value)  {conn[i] = value;};
        // GETTERS
        int  getConn (int index)         {return conn[index];};
};

/***************************************************************************************************
MESH DATA STRUCTURE
****************************************************************************************************
This class is used to keep the pointers to node and element data structures. It is more
convinient to pass the pointer to this class rather than passin the pointers to both element and
node arrays during function calls. So it is just for simplification.
***************************************************************************************************/
class triMesh
{
    private:
        // VARIABLES
        int ne;           // total number of elements
        int nn;           // number of nodes
        triNode*    node; // pointer to partition, node level data structure
        triElement* elem; // pointer to element level data structure

        double * xyz;

        //METHODS
        void readMeshFiles(inputSettings*);
        void swapBytes(char*, int, int);    // Used during file read operations
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        triMesh(){};
        // DESTRUCTOR
        ~triMesh()
        {
            delete[] node;    // Recovers memory for node level data structure.
            delete[] elem;    // Recovers memory for element level data structure.
        };
        // GETTERS
        int         getNe()            {return ne;};
        int         getNn()            {return nn;};

        double*     getXyz()           {return xyz;};

        triNode*    getNode(int index) {return &node[index];};
        triElement* getElem(int index) {return &elem[index];};

        // INTERFACE METHOD
        void prepareMesh(inputSettings*);
};


#endif /* TRI_H_ */
