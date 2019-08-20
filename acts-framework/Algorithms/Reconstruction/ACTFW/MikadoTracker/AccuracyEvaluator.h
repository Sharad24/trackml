#ifndef AccuracyEvaluator_H
#define AccuracyEvaluator_H


class Tracker;
class RecoPassParameters;

class AccuracyEvaluator
{
 public:
 
  AccuracyEvaluator( Tracker *tracker ) :mTracker( tracker) {}
  //AccuracyEvaluator() = default;  
  ~AccuracyEvaluator() = default;  
  
  void SelectParticles( const RecoPassParameters &par, double phiMin, double phiMax );
  void Evaluate( bool doPrint );  

  Tracker *mTracker=0;  

  double mEfficiency=-1;
  double mPurity=-1;
  double mFakes=-1;
  double mTotalEfficiency=-1;
  double mTotalRemoved=-1;  
  double mTotalFakes=-1;
  double mHitsMissed = 0;
  double mHitsFake = 0;

  
   int    statRecN[12]  ={0,0,0,0,0,0,0,0,0,0,0};
   int    statFakeN[12] ={0,0,0,0,0,0,0,0,0,0,0};
   int    statShortN[12]={0,0,0,0,0,0,0,0,0,0,0};
   double statRecW[12]  ={0,0,0,0,0,0,0,0,0,0,0};
   double statRecWRemoved[12]={0,0,0,0,0,0,0,0,0,0,0};
   double statFakeW[12]={0,0,0,0,0,0,0,0,0,0,0};
   double statShortW[12]={0,0,0,0,0,0,0,0,0,0,0};
   double statShortSelW[12]={0,0,0,0,0,0,0,0,0,0,0};
   double statRecWMissed = 0;
   double statRecWFake = 0;
   double statRecWAllHits = 0;
   double statPurity = 0.;

};

#endif
