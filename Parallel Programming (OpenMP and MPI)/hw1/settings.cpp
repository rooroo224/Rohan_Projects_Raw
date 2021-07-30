#include "settings.h"

/***************************************************************************************************
void inputSettings::inputSettings()
private - Default constructor
***************************************************************************************************/
inputSettings::inputSettings(int largc, char** largv)
{
    // Default values for the input parameters
    argc = largc;   //Command line parameters are stored at private var
    argv = largv;   //Command line parameters are stored at private var
    title = "Default Title";
    minfFile = "minf";
    mxyzFile = "mxyz";
    mienFile = "mien";
}

/***************************************************************************************************
void inputSettings::readSettingsFile()
public - main interface function, calls all other, private, functions of the class
***************************************************************************************************/
void inputSettings::prepareSettings()
{
    readSettingsFile();
    printSettings();

    return;
}


/***************************************************************************************************
void inputSettings::readSettingsFile()
private -reads settings.in file
***************************************************************************************************/
void inputSettings::readSettingsFile()
{
    char * fileName;
    string lineString;
    string dummyString;
    char   dummyChar;

    ifstream inputFile;
    if (argc > 1)
    {
        fileName = argv[1];
    }
    else
    {
        fileName = "settings.in";
    }

    inputFile.open(fileName,ios::in);
    if (inputFile.is_open()==false)
    {
        cout << "Unable to open input file: " << fileName << " ! Aborting... " << endl;
        exit(0);
    }

    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "minf")
                iss >> minfFile;
            else if(dummyString == "mxyz")
                iss >> mxyzFile;
            else if(dummyString == "mien")
                iss >> mienFile;
            else
            {
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                inputFile.close();
                exit(0);
            }
        }
    }
    inputFile.close();

    return;
}


/***************************************************************************************************
void inputSettings::printSettings()
***************************************************************************************************/
void inputSettings::printSettings()
{
    cout.precision(4);
    cout << scientific;
    cout << endl << "==================== SETTINGS ====================" << endl;
    cout << "Title of the simualation              : " << title        << endl;
    cout << "Name of the minf file                 : " << minfFile    << endl;
    cout << "Name of the mxyz file                 : " << mxyzFile    << endl;
    cout << "Name of the mien file                 : " << mienFile    << endl;

    return;
}


