#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "constants.h"

/***************************************************************************************************
Input parameters are stored in this class. readSettingsFile() method does the job: It opens the
settings.in file, reads the parameters and assigns them to the private variables.
Why triMEsh is a friend? : triMesh method readMeshFiles needs to set ne and nn. no other funtions
need to change a private variable of inputSettings. Thus, previliges are given to triMesh.
***************************************************************************************************/

class inputSettings
{
    private:
        // VARIABLES
        int    argc;     // command line argument
        char** argv;     // command line argument
        string title;    // title of the document
        string minfFile; // information file name
        string mxyzFile; // node coordinates file name
        string mienFile; // connctivity file name
        // METHODS
        void readSettingsFile();
        void printSettings();

    protected:

    public:
        /// CONSTRUCTORS ///
        inputSettings();
        inputSettings(int, char**);
        /// GETTERS ///
        int    getArgc()     {return argc;};
        char** getArgv()     {return argv;};
        string getTitle()    {return title;};
        string getMinfFile() {return minfFile;};
        string getMxyzFile() {return mxyzFile;};
        string getMienFile() {return mienFile;};

        //INTERFACE METHOD
        void prepareSettings();

};



#endif /* SETTINGS_H_ */
