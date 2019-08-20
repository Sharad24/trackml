#ifndef TrackSelector_H
#define TrackSelector_H


struct Tracker;

class TrackSelector
{
 public:

 TrackSelector( Tracker *tracker ) : mTracker( tracker ) {}
  ~TrackSelector() = default;

  Tracker *mTracker = 0;  
  void SelectTracks();

};




#endif
