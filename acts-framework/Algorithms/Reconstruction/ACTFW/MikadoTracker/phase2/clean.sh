#!/bin/bash

sudo rm -rf __pycache__ build *.so
rm -f *~ geoLayer* cuts.txt
cd testSubmission; . ./clean.sh; cd ..

