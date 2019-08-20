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

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy as np
import pandas as pd
import os
import ctypes 

from trackerInterface import trackerInterfaceInit as _trackerInterfaceInit
from trackerInterface import trackerInterfaceExec as _trackerInterfaceExec

__authors__ = ['Sabrina Amrouche', 'David Rousseau', 'Moritz Kiehn', 'Ilija Vukotic']

class Model():
    def __init__(self, eps=0.008, min_samples=1, metric='euclidean',
                 algorithm='kd_tree', leaf_size=30, p=None, n_jobs=1, verbose=1):
        super().__init__()
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.p = p
        self.n_jobs = n_jobs
        self.verbose = verbose
        path = os.path.dirname(os.path.realpath(__file__))
        _trackerInterfaceInit( path  );

    def predict_one_event(self, event_id, event, cells=None):        
        event1 = event.copy()
        x = event1.x.values
        y = event1.y.values
        z = event1.z.values
        ids = event1.hit_id.values
        vol = event1.volume_id.values
        layer = event1.layer_id.values
        labels = np.zeros(event1.x.size, np.int32)

        _trackerInterfaceExec(x, y, z, ids, vol, layer, labels)

        sub = pd.DataFrame(data=np.column_stack((event.hit_id.values, labels)),
                           columns=["hit_id", "track_id"]).astype(int)
        sub['event_id'] = event_id
        return sub
