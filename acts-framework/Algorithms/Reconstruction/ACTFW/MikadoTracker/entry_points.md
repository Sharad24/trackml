** To run the code standalone: ** 

1. specify the path to the data and N of events to process in reconstruction.cxx file. 
By default, it is "./data"

One can additionally set a flag in reconstruction.cxx for analysing the reconstruction efficiency with the truth data:
bool analyseTruth = true;

2. compile the code via:
make

3. run the reconstruction by typing:
./reco

the reconstruction output will be written to mysubmission.csv file


** To run the code with the official docker container: **

- To compile and run the code with the official competition docker container:

cd phase2

. compile.sh

. run.sh

- A .zip file for the official submission was created via:

cd phase2

. compile.sh

. pack.sh

- To test the submission.zip:

cd phase2

. compile.sh

. pack.sh

cd testSubmission

. run.sh

