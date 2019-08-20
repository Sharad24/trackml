#include "TrackModelPhysical.h"
#include <math.h>
#include <iostream>
#include "Tracker.h"


using namespace std;

int TrackModelPhysical::createXY(double xf, double yf, double zf,
				 double xm, double ym, double zm,
				 double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;//Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xm-= xf;
  ym-= yf;
  zm-= zf;

  xl-= xf;
  yl-= yf;
  zl-= zf;

  double L = sqrt( xl*xl + yl*yl );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv; 
  ey = yl*Linv; 


  double um =  xm*ex + ym*ey;
  double vm = -xm*ey + ym*ex;

  //SG!! if( um<0.1*L || um>0.9*L ) return -2;
  //if( um<0.1 || um>L-0.1 ) return -2;

  q = 1;
  pv = 0.5*L;

  if( vm<0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( um*(L-um) - vm*vm );
  double vmAbs = fabs(vm);
  if( tmp < 0.02*fabs(pv*vm) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vmAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vmAbs;

  int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  double sinA = 0.5*L*ptInv;
  if( sinA>.999 ) return -2;
  pz = 0.5*zl / asin(sinA);

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;

  return ret;
}


int TrackModelPhysical::createXYMiddle(double xf, double yf, double zf,
				       double xm, double ym, double zm,
				       double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xm;
  y0 = ym;
  z0 = zm;

  xf-= x0;
  yf-= y0;
  zf-= z0;

  xl-= x0;
  yl-= y0;
  zl-= z0;

  double L = sqrt( xl*xl + yl*yl );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv; 
  ey = yl*Linv; 

  double uf =  xf*ex + yf*ey;
  double vf = -xf*ey + yf*ex;

  q = 1;
  pv = 0.5*L;

  if( vf > 0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( uf*(-L+uf) + vf*vf );
  double vfAbs = fabs(vf);
  if( tmp < 0.02*fabs(pv*vf) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;

   int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  double sinA = 0.5*L*ptInv;
  if( sinA>.999 ) return -2;
  pz = 0.5*zl / asin(sinA);

  pu*=fabs(Bz);
  pv*=fabs(Bz);
  pz*=fabs(Bz);
  pt*=fabs(Bz);
  if( Bz < 0 ) q = -q;
  ptInv = 1./pt;

 return ret;
}


int TrackModelPhysical::createXYLast(double xf, double yf, double zf,
				     double xm, double ym, double zm,
				     double xl, double yl, double zl,
				     double Bz  )
{  
  const double radMax = 123./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  xf-= xl;
  yf-= yl;
  zf-= zl;

  xm-= xl;
  ym-= yl;
  zm-= zl;

  double L = sqrt( xm*xm + ym*ym );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  double ex1 = -xm*Linv; 
  double ey1 = -ym*Linv; 

  double uf =  xf*ex1 + yf*ey1;
  double vf = -xf*ey1 + yf*ex1;

  double tmp = 0.5* ( uf*(L+uf) + vf*vf );
  double vfAbs = fabs(vf);
  if( vfAbs < 1.e-4 ) return -2;//

  pv = -0.5*L;
  if( tmp < 0.02*fabs(pv*vf) ) return -3; // incl. angle > 88.8 grad
  
  x0 = xl;
  y0 = yl;
  z0 = zl;
  ex = ex1;
  ey = ey1;

  q = 1;

  if( vf > 0. ){
    pv = -pv;
    q = -1;    
  }
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;
  
  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  double sinA = 0.5*L*ptInv;
  if( sinA>.999 ) return -2;
  pz = -0.5*zm / asin(sinA);

  pu*=fabs(Bz);
  pv*=fabs(Bz);
  pz*=fabs(Bz);
  if( Bz < 0 ) q = -q;
  pt*=fabs(Bz);

  ptInv = 1./pt;

  return 0;
}


int TrackModelPhysical::createXYLastFix( double xm, double ym, double zm,
					 double xl, double yl, double zl,
					 double Bz, double fixPt, double fixQ  )
{
  //const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  pt = fixPt;
  q = fixQ;
  x0 = xl;
  y0 = yl;
  z0 = zl;
  
  double qBz = q*Bz;
  double rBz2 = pt*pt;

  xm-= xl;
  ym-= yl;
  zm-= zl;

  double L2 = xm*xm + ym*ym;
  if( L2 <= .1 ) return -1; 

  double t2 = rBz2/L2 - 0.25*Bz*Bz;
  if( t2 <= .1*Bz*Bz ) return -2; 

  double t = sqrt( t2 );
  
  double px = - (0.5*ym*qBz + xm*t);
  double py =   (0.5*xm*qBz - ym*t);

  double L = sqrt( L2 );
  double Linv = 1./L;

  ex = -xm*Linv;
  ey = -ym*Linv;

  pu =  px*ex + py*ey;
  pv = -px*ey + py*ex;
 
  ptInv = 1./pt;
  double sinA = 0.5*L*fabs(Bz)*ptInv;
  if( sinA >.999 ) return -3;
  pz = -0.5*zm / asin(sinA)*fabs(Bz);

  return 0;
}
 
int TrackModelPhysical::createXYFirst(double xf, double yf, double zf,
				      double xm, double ym, double zm,
				      double xl, double yl, double zl, double Bz  )
{
  const double radMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xl-= x0;
  yl-= y0;
  zl-= z0;

  xm-= x0;
  ym-= y0;
  zm-= z0;

  double L = sqrt( xm*xm + ym*ym );
  
  //if( L < 1. ) return -1; // 1cm cut on length
  if( L < .1 ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xm*Linv; 
  ey = ym*Linv; 

  double ul =  xl*ex + yl*ey;
  double vl = -xl*ey + yl*ex;

  q = 1;
  pv = 0.5*L;

  if( vl > 0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( ul*(-L+ul) + vl*vl );
  double vfAbs = fabs(vl);
  if( tmp < 0.02*fabs(pv*vl) ) return -3; // incl. angle > 88.8 grad
  
  if( tmp > radMax*vfAbs ){
    pu = radMax; // pu > 100 GeV -> pu = 100GeV 
    pv = 0;
  } else pu = tmp/vfAbs;

  pt = sqrt( pu*pu + pv*pv );
  //if( pt<.01 ) ret = -2;

  ptInv = 1./pt;
  double sinA = 0.5*L*ptInv;
  if( sinA>.999 ) return -2;
  pz = 0.5*zm / asin(sinA);

  pu*=fabs(Bz);
  pv*=fabs(Bz);
  pz*=fabs(Bz);
  pt*=fabs(Bz);
  if( Bz < 0 ) q = -q;
  ptInv = 1./pt;
 return 0;
}
 



int TrackModelPhysical::createZ(double xf, double yf, double zf,
				double xm, double ym, double zm,
				double xl, double yl, double zl, double Bz  )
{
  return createXYFirst(xf,yf,zf,xm,ym,zm,xl,yl,zl,Bz);


  const double puMax = 100./Geo::OriginBzkG;///Bz; // 100 GeV maximum for pu

  x0 = xf;
  y0 = yf;
  z0 = zf;

  xm-= xf;
  ym-= yf;
  zm-= zf;

  xl-= xf;
  yl-= yf;
  zl-= zf;

  double L = sqrt( xl*xl + yl*yl );
  
  if( L < 1. ) return -1; // 1cm cut on length

  double Linv = 1./L;

  ex = xl*Linv;
  ey = yl*Linv;

  double um =  xm*ex + ym*ey;
  double vm = -xm*ey + ym*ex;

  if( um<0.1*L || um>0.9*L ) return -2; //SG!! remove the cut later

  q = 1;
  pv = 0.5*L;

  if( vm<0. ){
    pv = -pv;
    q = -1;    
  }
      
  double tmp = 0.5* ( um*(L-um) - vm*vm );
  double vmAbs = fabs(vm);
  if( tmp < 0.02*fabs(pv*vm) ) return -3; // incl. angle > 88.8 grad  
  if( tmp > puMax*vmAbs ) pu = puMax; // pu > 100 GeV -> pu = 100GeV 
  else pu = tmp/vmAbs;

  int ret = 0;

  pt = sqrt( pu*pu + pv*pv );
  ptInv = 1./pt;
  //if( pt<.01 ) ret = -2;

  // calculate ptInv with respect to zm

  double cl = 0.5*L;
  double cm = 0.5*sqrt( um*um + vm*vm );
  for(int iter=0; iter<0; iter++){
    // we solve equation zm*asin(cl*ptInv) = zl*asin(cm*ptInv)
    // new ptInv = ptInv + d.
    double al = cl*ptInv;
    double am = cm*ptInv;
    // linearisation at d=0:
    // zm*( asin(al) + al/sqrt(1.-al*al)*d ) = zl*( asin(am) + am/sqrt(1.-am*am)*d )
    // d * [ zm*al/sqrt(1.-al*al) - zl*am/sqrt(1.-am*am) ] = zl*asin(am) - zm*asin(al)
    double bl = sqrt(1.-al*al);
    double bm = sqrt(1.-am*am);
    // d * [ (zm*al*bm - zl*am*bl)/(am*bm) ] = zl*asin(am) - zm*asin(al)
    double d = am*bm*( zl*asin(am) - zm*asin(al) ) / (zm*al*bm - zl*am*bl);
    ptInv+=d;
  }

  pt = 1./ptInv;
  pu = sqrt(pt*pt-pv*pv);
  pz = 0.5*zl/asin(0.5*L*ptInv);  

  pu*=Bz;
  pv*=Bz;
  pz*=Bz;
  pt*=Bz;
  ptInv = 1./pt;

  return ret;
}

int TrackModelPhysical::getDistanceAtXY( double x, double y, double z, double &duz, double &dv, double Bz, double *DS ) const
{
  dv = 0;
  duz = 0;
  x-=x0;
  y-=y0;
  z-=z0;
  
 
  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;
  
  double epv = pv - q*Bz*u;
  double epu2 = pt*pt - epv*epv;
  //SG!!! if( epu2<1 ) return -1;
  
  double epu = sqrt(epu2);
  double tv = ( pv + epv )/( pu + epu );

  double vTrack = u*tv;
  
  double chord = u*sqrt(1.+tv*tv); // chord to the extrapolated point == sqrt(u^2+vTrack^2)*sign(u)
  double sa = 0.5*chord*q*Bz*ptInv; //  sin( half of the rotation angle ) ==  (chord/2) / radius

  double dS;

  if( fabs(sa)>1.e-2 && fabs(sa)<.999){
    dS = q*pt/Bz*2.*asin( sa ); // path in XY
  } else {
    dS = q*pt/Bz*2.*sa*(1. + (1./6.)*sa*sa); // path in XY
  }

  double dLp = ptInv*dS; // path in XYZ / p == path in XY / pt
  double zTrack = pz * dLp;
  
  dv = v-vTrack;
  duz = z-zTrack;
  if( DS ) *DS=fabs(u);


  return 0;
}


int TrackModelPhysical::getDistanceAtXY( double alpha, double x, double y, double z, double &duz, double &dv, double Bz ) const
{
  dv = 0;
  duz = 0;
  x-=x0;
  y-=y0;
  z-=z0;
 
  // get px, py:  rotate pu,pv to -(ex,ey)
  
  double px =  pu*ex - pv*ey;
  double py = +pu*ey + pv*ex;

  double sA=sin(alpha);
  double cA=cos(alpha);

  // rotate x axis to alpha

  double u =  x*cA + y*sA;
  double v = -x*sA + y*cA;
  double pur =  px*cA + py*sA;
  double pvr = -px*sA + py*cA;
  

  double epv = pvr - q*Bz*u;
  double epu2 = pt*pt - epv*epv;
  //SG!!! if( epu2<1 ) return -1;
  
  double epu = sqrt(epu2);
  double tv = ( pvr + epv )/( pur + epu );

  double vTrack = u*tv;
  
  double chord = u*sqrt(1.+tv*tv); // chord to the extrapolated point == sqrt(u^2+vTrack^2)*sign(u)
  double sa = 0.5*chord*q*Bz*ptInv; //  sin( half of the rotation angle ) ==  (chord/2) / radius
  if( fabs(sa)>.999 ) return -1;
  double dS = q*pt/Bz*2.*asin( sa ); // path in XY
  double dLp = ptInv*dS; // path in XYZ / p == path in XY / pt
  double zTrack = pz * dLp;

  dv = v-vTrack;
  duz = z-zTrack;
  return 0;
}



int TrackModelPhysical::getDistanceAtZ( double x, double y, double z, double &duz, double &dv, double Bz, double *DS ) const
{
  duz = 0;
  dv = 0;
  x-=x0;
  y-=y0;
  z-=z0;
 
  double tz = z/pz;
  double a = tz*Bz; // = s/pt
  double f1, f2;
  if( fabs(a)>1.e-2 ){
    f1 = sin(a)/Bz;
    f2 = q*(1.-cos(a))/Bz;
  } else {
    f1 = tz*(1. - a*a*(1./6.)) ;
    f2 = q*tz*0.5*a;  
  }

  double uTrack = (f1*pu + f2*pv);
  double vTrack = (f1*pv - f2*pu);      
  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;
  duz = u - uTrack;
  dv  = v - vTrack;
  if( DS ){
    *DS=fabs(uTrack);
    //*DS=fabs(z);
  }
 
  return 0;
}



int TrackModelPhysical::getPhiZatR( double r, double &phi, double &z, double Bz ) const
{

  if( pt < 20 ){
    // TODO: simplify
    
    // xc,yc center of circle
    
    double px = (pu*ex - pv*ey)/Bz;
    double py = (pu*ey + pv*ex)/Bz;

    double xc =  x0 + py*q;
    double yc =  y0 - px*q;
    
    double rc = sqrt(xc*xc+yc*yc);
    //cout<<"circle xc "<<xc<<" yc "<<yc<<" rc "<<rc<<endl;  
    //cout<<"mx 1"<<endl;
    if( rc<1 ) return  -1; // center of helix is too close to the origin
    
    double rci = 1./rc;
    double r1 = r*rci;
    double pt1 = (pt/Bz)*rci;
    double a = 0.5*( 1.  + r1*r1 - pt1*pt1 );
    
    double b = r1*r1-a*a;
    //cout<<"mx 2: b = "<<b<<endl;
    if( b<1.e-8 ){
      //Print();
      //cout<<"xc = "<<xc<<" yc = "<<yc<<endl;
      return -2;
    }
    b = sqrt(b);
    
    double x = xc*a - yc*b;
    double y = yc*a + xc*b;
    double chord = (x-x0)*(x-x0)+(y-y0)*(y-y0);

    double x2 = xc*a + yc*b;
    double y2 = yc*a - xc*b;
    double chord2 = (x2-x0)*(x2-x0)+(y2-y0)*(y2-y0);
  
    if( chord > chord2 ){
      x = x2;
      y = y2;
      chord = chord2;
    }
  
    phi = atan2(y,x);
  

    double u = (x-x0)*ex + (y-y0)*ey;
    //double v =-(x-x0)*ey + (y-y0)*ex;  
 
    // get z at u:
  
    chord = sqrt(chord); // chord to the extrapolated point == sqrt(u^2+vTrack^2)*sign(u)

    if( u<0 ) chord = -chord;
    double sa = 0.5*chord*q*(ptInv*Bz); //  sin( half of the rotation angle ) ==  (chord/2) / radius
    // double dS = q*(pt/Bz)*2.*asin( sa ); // path in XY
    // double dLp = (ptInv*Bz)*dS; // path in XYZ / p == path in XY / pt
    if( fabs(sa)>.999 ) return -2;
    double dLp = q*2.*asin( sa ); // path in XYZ / p == path in XY / pt
    z = z0 + (pz/Bz) * dLp;
    //cout<<"mx 3"<<endl;

  } else {
    //cout<<"mx 4"<<endl;

    double a = x0*ey-y0*ex;
    double d = r*r - a*a;
    if( d<1.e-2 ) d = 0.;
    d = sqrt(d);
    double b = x0*ex+y0*ey;
    double t = -b+d;
    double t1 = -b-d;    
    if( fabs(t1)<fabs(t) ) t = t1;
    double x = x0 + t*ex;
    double y = y0 + t*ey;
    phi = atan2(y,x);
    z = z0 +  t*pz*ptInv;    
    //cout<<"mx 5"<<endl;
  }

  return 0;
}



int TrackModelPhysical::getPhiRatZ( double z, double &phi, double &r, double Bz ) const
{
  phi = 0.;
  r = 0.;

  if( fabs(pz)<1.e-4 ) return -1; 

  z-=z0;     
  double qBz = q*Bz;
  double tz = z/pz;
  
  // rotation angle
  double a = -tz*qBz; // = s/pt

  double f1,f2;
  if( fabs(a)>1.e-2 ){
    f1 = -sin(a)/qBz;
    f2 = (1.-cos(a))/qBz;      
  } else {
    f1 = tz*(1. - a*a*(1./6.)) ;
    f2 = -tz*0.5*a;  
  }
  double u = (f1*pu + f2*pv);
  double v = (f1*pv - f2*pu); 
  double x = x0 + u*ex - v*ey;
  double y = y0 + u*ey + v*ex;      
   
  r = sqrt(x*x+y*y);
  phi = atan2(y,x);

  // check max rotation
  {  
     // xc,yc center of circle    
    double px = (pu*ex - pv*ey);
    double py = (pu*ey + pv*ex);
 
    double dr = ( x0*px + y0*py  )*tz;

    double xc =  x0 + py / qBz;
    double yc =  y0 - px / qBz;
    
    double rc2 = xc*xc+yc*yc; 
    if( rc2<1 ) return  -1; // center of helix is too close to the origin 

 
    // rot. angle of the outermost point of the circle

    double a0 = atan2( y0 - yc, x0 - xc );
    double amax0 = 0;
    if( dr>=0 ){
      amax0 = atan2(  yc,  xc );
    } else {
      amax0 = atan2(  -yc,  -xc );
    }
    double amax = fabs(amax0 - a0);
    if( amax >  2*M_PI ) amax -= 2*M_PI;    
    if( amax >  M_PI ) amax = 2*M_PI - amax;    
    /*
    cout<<"a0 = "<<a0/M_PI*180.<<" amax0 = "<<amax0/M_PI*180.<<endl;
    cout<<"a = "<<a/M_PI*180.<<" amax = "<<amax/M_PI*180.<<endl;
    cout<<" extr z = "<<z0+z<<" tz = "<<tz<<endl;
    cout<<"xc "<<xc<<" yc "<<yc<<endl;
    cout<<" q "<<q<<" Bz "<<Bz<<endl;
    Print();
    */        
    if( fabs(a) >= amax ){
      //cout<<"skip layer!!!"<<endl;
      return -1234;
    }
  }    

  return 0;
}


void TrackModelPhysical::Print() const
{
  cout<<" x "<< x0<<" y "<<y0<<" z "<<z0;
  cout<<" ex "<<ex<<" ey "<<ey<<" pt "<<pt<<" pz "<<pz<<" q "<<q<<endl;
}

int TrackModelPhysical::estmateBzkG( const HitMC &mc, double x, double y, double &BzkG )
{  
  if( mc.q!=1. && mc.q!=-1. ){
    cout<<"estimate Bz: wrong charge!!!"<<endl;
    exit(1);
  }
  double ex = mc.px/mc.pt; 
  double ey = mc.py/mc.pt;
  
  x-=mc.x;
  y-=mc.y;

  double u =  x*ex + y*ey;
  double v = -x*ey + y*ex;
  double l2 = u*u + v*v;

  if( l2<.3*.3 ) return -1; // require at least 3 mm distance in XY

  BzkG = (-2*mc.q*mc.pt*v) / l2 / Geo::CLight;
  return 0;
}


