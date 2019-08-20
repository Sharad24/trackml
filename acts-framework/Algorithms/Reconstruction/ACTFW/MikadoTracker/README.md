Hi!  

Below you can find an outline of how to reproduce my solution for the TrackML competition.  
If you run into any trouble with the setup/code or have any questions please contact me at sergey.gorbuov.32@gmail.com  

The algorithm is described here:  [doc/MikadoTrackerForTrackML.pdf](https://gitlab.com/sgorbuno/MikadoTracker/blob/master/doc/MikadoTrackerForTrackML.pdf)


#ARCHIVE CONTENTS  
Makefile            : file which steers the compilation  
*.cxx *.h           : reconstruction code  
analysis/*.cxx/*.h  : additional code which is used to analyse the data  
cuts.txt,
geoLayerField.txt,  
geoLayerSizes.txt   : configuration files with detector geometry and magnetic field approximation  


#HARDWARE: (The following specs were used to create the original solution)  
Mac OsX 10.13.6  
2,6 GHz Intel Core i5  
8 GB 1600 MHz DDR3  


#SOFTWARE:  

(optional) ROOT 6.14/00 or any other version of ROOT   
"root-config" script should be accesible in the command line   
if it is not present, the code will be compiled without ROOT

Web-site of the package: https://root.cern.ch/   
Installation instructions can be found here: https://root.cern.ch/downloading-root  

The package is only used to analyse the data. 


#DATA SETUP (assumes the [Kaggle API](https://github.com/Kaggle/kaggle-api) is installed)  
 below are the shell commands used in each step, as run from the top level directory  

cd data/  
kaggle competitions download -c trackml-particle-identification -f train_sample.zip  
kaggle competitions download -c trackml-particle-identification -f test.zip  
unzip train_sample.zip  
unzip test.zip  

One event is already present in the "./data/" folder.

#DATA PROCESSING FOR DEVELOPMENT 

1. specify the path to the data and N of events to process in reconstruction.cxx file.   
By default, it in "./data/"  

One can additionally set a flag for analysing the reconstruction efficiency with the truth data.
In the reconstruction.cxx:  

bool analyseTruth = true;  

2. compile the code via:  
make  

3. run the reconstruction by typing:  
./reco  

the reconstruction output will be written to mysubmission.csv file  

#DATA PROCESSING WITH THE OFFICIAL DOCKER CONTAINER

To compile and run the code with the official competition docker container:

cd phase2
. compile.sh
. run.sh

A .zip file for the official submission was created via:

cd phase2
. compile.sh
. pack.sh

To test the submission.zip:
cd phase2
. compile.sh
. pack.sh
cd testSubmission
. run.sh

#MODEL BUILD:  

The algorithm is combinatorial, there is no usual model build step.   
But it has some parameters, which are extracted from the test samples and stored in    
cuts.txt, geoLayerField.txt,  geoLayerModules.txt,  geoLayerSizes.txt files.   

In order to rebuid the geometrical parameters, one should uncomment   
  //tracker.analyzeGeometry(0);     
  //continue;  
and  
  //tracker.analyzeGeometry(1);  
lines in the reconstruction.cxx  

Reconstruction parameters are stored in cuts.txt. 
In order to tune the parameters  one should specify a reconsrtruction pass in reconstruction.cxx :

 tracker.mRecoPassForLearning = /*pass number*/;

 The new parameters will be saved in cutsNew.txt file, which can then replace cuts.txt