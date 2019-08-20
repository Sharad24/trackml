#ifndef EventReader_H
#define EventReader_H

#include "Geo.h"
#include "DataStructures.h"

class Tracker;

class EventReader
{
 public:

 EventReader( Tracker *tracker ) : mTracker( tracker ) {}
  ~EventReader() = default;

  int readEvent( const char *directory, int event, bool loadMC );
  int readEvent( float *x, float *y, float *z, int *id, int *volumes, int *layers, int nHits );  
  int readMCEvent( const char *directory, int event );
  bool checkLayers( const int *baseLayers, int nBaseLayers, Particle &p, double &fitpt, double &phi );

  Tracker *mTracker = 0;  
  int layerNHits[Geo::NLayers];
};




#endif
