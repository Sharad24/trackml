// Neural Network based tracker
// M.Kunze, Heidelberg University, 2018

#include "Tracker.h"
//#include "Python.h"
#include <iostream>
#include <string>

Tracker *tracker = NULL;

extern "C" {
    void processSetup(int filenum=0,const char *datapath=DATAPATH,const char *workpath=WORKPATH)
    {
        std::cout << "processSetup " << filenum << " " << datapath << " " << workpath << std::endl;
        if (tracker != NULL) delete tracker;
        tracker = new Tracker(filenum,datapath,workpath);
        //Py_INCREF(tracker);
    }
    void processInitHits(int nhits,double *x,double *y,double *z,int* v,int* l,int* m)
    {
        std::cout << "processInitHits " << nhits << std::endl;
        tracker->initHits(nhits,x,y,z,v,l,m);
    }
    void processInitCells(int ncells,int *id,int *ch0,int *ch1,double *value)
    {
        std::cout << "processInitCells " << ncells << std::endl;
        tracker->initCells(ncells,id,ch0,ch1,value);
    }
    int *processFindTracks(int step,const char *path)
    {
        std::cout << "processFindTracks " << step << " " << path << std::endl;
        int *labels = tracker->findTracks(step, path);
        return labels;
    }
    void processReadTruth()
    {
        tracker->readTruth();
    }
    void processReadStarts()
    {
        tracker->readStarts();
    }
    void processReadHits()
    {
        tracker->readHits();
    }
    void processReadCells()
    {
        tracker->readCells();
    }
    void processReadBlacklist()
    {
        tracker->readBlacklist();
    }
    void processReadWhitelist()
    {
        tracker->readWhitelist();
    }
    void processSortTracks(void)
    {
        std::cout << "processSortTracks " << std::endl;
        tracker->sortTracks();
    }
    void processFinish(void)
    {
        std::cout << "processFinish " << std::endl;
        //Py_DECREF(tracker);
        delete tracker;
        tracker = NULL;
    }

}

