#ifndef TRACKMODELPHYSICAL_H
#define TRACKMODELPHYSICAL_H

#include "Geo.h"
#include <limits>

struct HitMC;

class TrackModelPhysical
{
 public:
  //TrackModelPhysical(){ init();}
  
  void initInf(){
    init( std::numeric_limits<double>::infinity() );
  }

  void init( double v = 0.){    
    x0=v; // position of the parameterization point
    y0=v;
    z0=v;
    ex=v; // directions of u axis
    ey=v;
    pu=v; // momentum
    pv=v;
    pz=v; 
    q=v;
    pt=v;  
    ptInv=v;
  }

  double x0; // position of the parameterization point
  double y0;
  double z0;
  double ex; // directions of u axis
  double ey;
  double pu; // momentum
  double pv;
  double pz; 
  double  q;
  double pt;  
  double ptInv;

  //

  double r0() const { return sqrt(x0*x0+y0*y0); }

  int createXY(double xf, double yf, double zf,
	       double xm, double ym, double zm,
	       double xl, double yl, double zl, double Bz  );

  int createXYFirst(double xf, double yf, double zf,
		    double xm, double ym, double zm,
		    double xl, double yl, double zl, double Bz  );

  int createXYMiddle(double xf, double yf, double zf,
		     double xm, double ym, double zm,
		     double xl, double yl, double zl, double Bz  );

  int createXYLast(double xf, double yf, double zf,
		   double xm, double ym, double zm,
		   double xl, double yl, double zl, double Bz );

  int createXYLastFix( double xm, double ym, double zm,
		    double xl, double yl, double zl,
		    double Bz, double fixPt, double fixQ  );


  int createXYLine(double xf, double yf, double zf,		   
		   double xl, double yl, double zl  );

  int createZ( double xf, double yf, double zf,
	       double xm, double ym, double zm,
	       double xl, double yl, double zl, double Bz  );

  int getDistanceAtXY( double x, double y, double z, double &duz, double &dv, double Bz, double *DS=0 ) const;
  int getDistanceAtXY( double alpha, double x, double y, double z, double &duz, double &dv, double Bz ) const;

  int getDistanceAtZ( double x, double y, double z, double &duz, double &dv, double Bz, double *DS=0) const;
  
  int getDistanceAt( bool isRadialLayer, double x, double y, double z, double &duz, double &dv, double Bz ) const;
 
  double getPtkG() const { return pt;}//*Geo::OriginBzkG; }
  
  int getPhiZatR( double r, double &phi, double &z, double Bz ) const;
  int getPhiRatZ( double z, double &phi, double &r, double Bz ) const;
 
  int getPhiT( bool isRadialLayer, double z, double r, double &phi, double &t, double Bz ) const;

  double getP(){ return sqrt(pt*pt+pz*pz); }


  void Print() const;

  static int estmateBzkG( const HitMC &mc, double x, double y, double &BzkG );

};
 
inline int TrackModelPhysical::getPhiT( bool isRadialLayer, double z, double r, double &phi, double &t, double Bz ) const
{
  if( isRadialLayer ) return getPhiZatR( r, phi, t, Bz );
  else return getPhiRatZ( z, phi, t, Bz );
}

inline int TrackModelPhysical::getDistanceAt( bool isRadialLayer, double x, double y, double z, double &duz, double &dv, double Bz ) const
{
  if( isRadialLayer ) return getDistanceAtXY( x, y, z, duz, dv, Bz );
  else return getDistanceAtZ( x, y, z, duz, dv, Bz);
}



inline int TrackModelPhysical::createXYLine(double xf, double yf, double zf,				  
					    double xl, double yl, double zl )
{
  const double ptMax = 100.; // 100 GeV maximum for pu
  const double ptMaxInv = 1./ ptMax;
  
  x0 = xl;
  y0 = yl;
  z0 = zl;

  xl-= xf;
  yl-= yf;
  zl-= zf;
  
  double et = 1./sqrt(xl*xl+yl*yl);
  ex = xl*et;
  ey = yl*et;
  
  q = 1;
  pu = ptMax; // pu > 100 GeV -> pu = 100GeV 
  pv = 0;  
  pt = ptMax;
  ptInv = ptMaxInv;
  pz = zl*et*ptMax;
  return 0;
}
 
#endif
