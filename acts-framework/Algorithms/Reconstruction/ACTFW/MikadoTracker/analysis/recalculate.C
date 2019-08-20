#include <iomanip>
using namespace std;
#include "../FitCorrections.h"
 

int recalculate()
{


  Correction h0[9]=
    {
      { 0, 1 },
      { 0, 1 },
      { 0, 1 },

      { 0, 1 },
      { 0, 1 },
      { 0, 1 },

      { 0, 1 },
      { 0, 1 },
      { 0, 1 }
    };
  
 
  Correction hUZ[9]=
    {
      { 0, 0.9994},
      { 0, 1},
      { 0, 1 },

      { 0, 1},
      { 0, 1},
      { 0, 1 },

      { 0, 1},
      { 0, 1},
      { 0, 1 }
    };

  Correction hV[9]=
    {
      { 0, 1},
      { 0, 1},
      { 0, 1 },

      { 0, 1},
      { 0, 1},
      { 0, 1 },

      { 0, 1 },
      { 0,  1},
      { 0, 1 }
    };

   Correction hPickUZ[9]=
    {
      { 0, 1 },
      {  0, 1},
      { 0, 1 },

      { 0, 1 },
      { 0, 1  },
      { 0, 1 },

      { 0, 1 },
      { 0,  1},
      { 0, 1 }
    };

  Correction hPickV[9]=
    {
      { 0, 1 },
      { 0, 1},
      { 0, 1 },

      { 0, 1 },
      { 0.0004273, 1.002},
      { 0, 1 },

      { 0, 1},
      { 0, 1 },
      { 0, 1 }
    };

  updateCorrection("corrMsUZ", corrMsUZ, h0);
  updateCorrection("corrMsV", corrMsUZ, h0);

  updateCorrection("corrFitUZ", corrFitUZ, hUZ);  
  updateCorrection("corrFitV", corrFitV, hV);
  
  updateCorrection("corrPickUZ", corrPickUZ, hPickUZ);
  updateCorrection("corrPickV", corrPickV, hPickV);
  
 
  return 0;
}
