#include "DataStructures.h"
#include "Geo.h"

using namespace std;

void Hit::Print() const
{
  cout<<"hit: vol "<<Geo::getVolume(layerID)<<" l "<<Geo::getVolumeLayer(layerID)<<" x "<<x()<<" y "<<y()<<" z "<<z()<<" r "<<r<<" phi "<<phi<<endl; 
}

void HitMC::Print() const
{  
  cout<<"hitmc: "<<partID<< "  x "<<x<<" y "<<y<<" z "<<z<<" p "<<p<<" pt "<<pt<<" pz "<<pz<<" p "<<p<<endl; 
}

void  Particle::Print() const
{
  cout<<"particle: count id "<<countID<<" x "<<x<<" y "<<y<<" z "<<z<<" pt "<<pt<<" pz "<<pz<<" p "<<p<<endl;
};
