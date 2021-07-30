/***************************************************************************************************
Name        : hw2.cpp
Description : This code is designed for educational purposes.
              - Implement one-sided MPI communication
                - MPI_Get
                - MPI_Accumulate
***************************************************************************************************/

#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

void localizeArray(double* localArray, double* globalArray,int lower, int upper, int nnc, int mnc, int mype, int npes)
{
    // Add code here
  // Window is created within COMM_WORLD communicator
  // so that the data in the winow can be accessed by all the other processes
    MPI_Win win;
    MPI_Win_create(localArray, nnc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    // synchronysation point
    MPI_Win_fence(0,win);
    // MPI_Get is used to combine the data from local arrays into individual global array
    MPI_Get((globalArray + mype*mnc),nnc,MPI_DOUBLE,mype,0,nnc,MPI_DOUBLE,win);
    for(int i=0;i<npes;i++)
      {
	//implement  corner case, i.e for the process in which the size of array is smaller than mnc
	if(i!=mype)
	  {
	    if( nnc==mnc)
	  MPI_Get((globalArray + i*mnc),nnc,MPI_DOUBLE,i,0,nnc,MPI_DOUBLE,win);
	    if(nnc!=mnc)
	  MPI_Get((globalArray + i*mnc),mnc,MPI_DOUBLE,i,0,mnc,MPI_DOUBLE,win);
	  }
      }
    MPI_Win_fence(0,win);
    // the MPI window created is freed from memory space
    MPI_Win_free(&win);
}

void accumulateArray(double* localArray,double* sumArray, int nnc,int mnc,int mype,int npes )
{
    // Add code here
  MPI_Win win;
  // Window is created of all the localArrays in the communicator
  MPI_Win_create(localArray, nnc*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
  MPI_Win_fence(0,win);
  for (int j=0;j<npes;j++)
    {
      if(j!=mype)
	{
	  //Accumulate operation is performed to obtain the sum of local array of the
	  //particular process with all other local Arrays
	  MPI_Accumulate(localArray,nnc,MPI_DOUBLE,j,0,nnc,MPI_DOUBLE,MPI_SUM,win);
	}
    }
   MPI_Win_fence(0,win);
   // the obtained sums are passes into sumAarry for each of the process
   MPI_Get(sumArray,nnc,MPI_DOUBLE,mype,0,nnc,MPI_DOUBLE,win);
   MPI_Win_fence(0,win);
   MPI_Win_free(&win);
}

double checkSum(double* arr, int val)
{
    // Add code here
  double pass=0;
  // return of sum of all elements in the array 
  for (int k=0; k<val; k++)
    pass=pass+arr[k];
  return pass;
}

void displayResults(double* localArray, double* globalArray, double* sumArray ,int nn, int nnc, int mnc, int mype)
{
    cout << "##################################################################" << endl;
    cout << endl;
    cout << "mype:" << mype << " nnc:" << nnc << " mnc:" << mnc << endl;
    cout << "mype: " << mype << " localArray:" << endl;
    cout << "[ ";
    for(int i=0; i<nnc; i++){
        cout << localArray[i] << " " ;
    }
    cout << " ]" << endl;

    cout << "mype: " << mype << " globalArray:" << endl;
    cout << "[ ";
    for(int i=0; i<nn; i++)
    {
        if (mype * mnc == i) cout << "-> ";
        if (mype * mnc + nnc == i) cout << "<- ";
        cout << globalArray[i] << " ";
    }
    cout << "]" << endl;

    cout << "mype: " << mype << " sumArray[ ";
    for(int i=0; i<nnc; i++)
    {
        cout << sumArray[i] << " ";
    }
    cout << "]" << endl;

    cout << endl;
    cout << "mype: " << mype << " checkSum globalArray: " << checkSum(/* Add code here */ globalArray, nn) << endl;
    cout << "mype: " << mype << " checkSum sumArray: " << checkSum(/* Add code here */sumArray,nnc) << endl;
    cout << endl;

}

int main(int argc, char **argv) {
/***************************************************************************************************
MAIN PROGRAM FLOW
1. Initialize the data
2. MPI_Get
3. MPI_Accumulate
4. Display results
***************************************************************************************************/

    int mype, npes;
    double starttime;
    int nn, nnc, mnc;
    double * localArray;
    double * globalArray;
    double * sumArray;
    //extra array is created to store the value of localArray temporatily to Display 
    // as localArray is modified in Accumulate function
    double * extralocal; 

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    // 1. Define length of the main vector
    nn=10;
    // 2. Determine nnc, mnc

    // nnc is calculated by taking upperbound of ratio of size 
    // of main vector and number of processors
    nnc=ceil((double)nn/npes);
    mnc = nnc;
    // this is to account the corner case, in which the size of  
    // last partition is smaller,so that nnc will has to updated
	if ((mype+1) * mnc > nn)
	{
	  nnc = nn - mype * mnc;
	  if (nnc < 0)
	    {
	      nnc = 0;
	    }
	}
    // 3. Create localArray, globalArray and sumArray
    // Dynamic memory allocation for all the three vectors
	localArray  = new double[nnc*sizeof(double)];
        sumArray    = new double[nnc*sizeof(double)];
        globalArray = new double[nn*sizeof(double)];
	extralocal  = new double[nnc*sizeof(double)];

    // 4. Initialize the data in the different array.
    // localArray has to be initialized with random values between 0 and 10
	srand(mype+1);
    //here the random values are generated using srand() & rand() within the range
    //these are stored in the local array
    int upper=10; int lower=0;
	for(int i=0;i<nnc;i++)
	{
	  localArray[i] = (upper-lower)*((double)rand() / (double)RAND_MAX) +lower;
	  extralocal[i] = localArray[i]; //copying the values for display

	}

    // 5. Transfer data from remote to local (MPI_Get)
	localizeArray(localArray,globalArray,lower,upper,nnc,mnc,mype,npes);

    // 6. Accumulate data from local to remote (MPI_Accumulate)
	accumulateArray(localArray,sumArray,nnc,mnc,mype,npes);

    // 7. Display the resulting arrays together with their checksum
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < npes; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (i == mype) {
	  displayResults(extralocal, globalArray,sumArray, nn, nnc, mnc, mype);
        }
    }
    
    MPI_Finalize();

    return 0;
}
