#!/usr/bin/env python
# coding: utf-8

########################################################################
# ======================  TrackML CHALLENGE MODEL  =====================
########################################################################
# Author: Isabelle Guyon, Victor Estrade
# Date: Apr 10, 2018

# ALL INFORMATION, SOFTWARE, DOCUMENTATION, AND DATA ARE PROVIDED "AS-IS".
# PARIS-SUD UNIVERSITY, THE ORGANIZERS OR CODE AUTHORS DISCLAIM
# ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE, AND THE
# WARRANTY OF NON-INFRIGEMENT OF ANY THIRD PARTY'S INTELLECTUAL PROPERTY RIGHTS.
# IN NO EVENT SHALL PARIS-SUD UNIVERSITY AND/OR OTHER ORGANIZERS BE LIABLE FOR ANY SPECIAL,
# INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF SOFTWARE, DOCUMENTS, MATERIALS,
# PUBLICATIONS, OR INFORMATION MADE AVAILABLE FOR THE CHALLENGE.
from distutils.core import setup
from distutils.core import Extension
import numpy

# Create one module for each Python/C module you want to build
module = Extension('trackerInterface', 
                   [
                       'trackerInterface.c', 
                       'trackerCXXInterface.cxx',
                       '/home/tracker/AccuracyEvaluator.cxx',
                       '/home/tracker/Cuts.cxx',
                       '/home/tracker/DataStructures.cxx',
                       '/home/tracker/Engine.cxx',
                       '/home/tracker/EventReader.cxx',
                       '/home/tracker/Geo.cxx',
                       '/home/tracker/Learning.cxx',                       
                       '/home/tracker/SearchLayer.cxx',
                       '/home/tracker/Tracker.cxx',
                       '/home/tracker/TrackModelPhysical.cxx',
                       '/home/tracker/TrackSelector.cxx'
                   ],
                   extra_compile_args=['-Ofast'],  # Here you can add extra compilation options
#                   libraries = ['test'],
#                   library_dirs = ['/home/code'],
                   include_dirs=[numpy.get_include(),'/home/tracker'],  # include numpy C libraries
                   )
# setup() will run the compilation commands
setup(
    name="trackerInterface",
    version="0.1",
    description="Interface to c++ tracking code.",
    ext_modules=[module],  # The module list, here there is only one.
)
