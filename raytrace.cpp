/* ClusterCaustics
 no recenter on particle files
 */

#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include "gridmap.h"
#include <math.h>
#include "particle_halo2.h"
#include "utilities.h"
#include "concave_hull.h"
//#include "gadget.hh"

using namespace std;

int main(int arg,char **argv){
  
  time_t t;
  time(&t);
  
  const std::string output_dir = "OutputTest3/";
  Utilities::LOGPARAMS logparams(output_dir + "logfile.txt");

  COSMOLOGY cosmo(CosmoParamSet::Planck18);
  
  cosmo.setOmega_matter(0.24,true);
  cosmo.sethubble(0.72);
  
  const int Npix =  1000;// 1D number of pixels in initial grid of rays
  logparams("Npix",Npix);
  const int Nsmooth = 5;// nearest neighbors smoothing scale, 5 is default
  logparams("Nsmooth",Nsmooth);

  // A spherical region of the simulation can be extracted
  // radius of extracted region in comoving pc/h, if <=0 the whole simulation box is used
  const double Rmax = 1.0e5; // 100 kpc/h
  logparams("Rmax",Rmax);
  // set position around which you want to extract particles
  // this must be in the units of the simulation output.  In this case
  // it is comoving pc/h.
  Point_3d<float> target(500689.656250 , 498223.187500, 494912.468750);
  
  //const bool cluster_on = true;  // toggle including cluster simulation
  const bool do_maps = true;    // map output maps or no

  // multiple source redshifts planes can be set
  //const std::vector<double> zss = {2.0,1.5,1.0,0.9,0.8,0.7,0.6,0.5};         // source redshift
  std::vector<double> zss = {2.5};         // source redshift
  logparams("zss",zss);

  //long seed = -11920;
  long seed = time(&t);
  
  // which projection of simulation
  // specifies which axis is the line of sight
  // 1 = x-axis, 2 = y-axis, 3 = z-axis
  const int projection = 3; // 1,2 or 3
  logparams("projection",projection);
  
  // Path and name of Gadget2 snapshot file 
  const std::string inputfilename = "DataFiles/snap_058";
  logparams("inputfilename",inputfilename);

  // prefix put on output file names
  const std::string fileprefix = "snap_058";
  logparams("fileprefix",fileprefix);

  if(!Utilities::IO::check_directory(output_dir)){
    std::cout << "Creating directory " << output_dir << std::endl;
    Utilities::IO::make_directories(output_dir);
  }

  // make the rotation angle for the projection
  Point_2d theta;
  
  if(projection>=1 && projection<=3){
    if(projection==1) {
      theta[0]=0;
      theta[1]=PI/2;
    }
    if(projection==2) {
      theta[0]=PI/2;
      theta[1]=PI/2;
    }
    if(projection==3) {
      theta[0]=0;
      theta[1]=0;
    }
  }
  
  // Here a vector of LensHaloParticles is created.  
  // In this case the Gadget2 snapshot contains 6 types of particles,
  // so the vector will contain 6 LensHaloParticles objects.  
  // The smoothing is done according to the density of each type of particle in 3D. 
  time_t start_time = time(nullptr);
  std::vector<LensHaloParticles<float> > halos = LensHaloParticles<float>::make(
      inputfilename,SimFileFormat::gadget2    
      ,1.0e5   /// masses scale, this times value in file is in solar masses/h
      ,1.0     /// length scale, this times value in file is in pc/h
      ,1.0e-13 /// inverse area for mass compensation
      ,cosmo   /// cosmology
      ,Nsmooth /// number of neighbours for adaptive smoothing
      ,8       /// buckets size in tree
      ,0.1     /// opening angle for tree
      ,Rmax  /// radius of extracted region in comoving kpc/h, if <=0 no excision is done
      ,target  /// center of extracted region in comoving kpc/h
      ,true);
 
  time_t end_time = time(nullptr);
  std::cout << "times to construct halos: " << difftime(end_time,start_time) << " seconds" << std::endl;
  logparams("time to construct halos",difftime(end_time,start_time));

  target *= 0;  // simulation is recentered on target

  // gets the original redshift that was in the sim file
  // This could be changed with halos[i].setRedshift();
  double zl = halos[0].getZlens();
  
  // angular size distance to the lens redshift
  double Dl = cosmo.angDist(zl);
  
  std::cout << "zs = " << zss[0] << std::endl;
  std::cout << "zl = " << zl << std::endl;
  std::cout << "Dls/Ds = " << cosmo.angDist(zl,zss[0])/cosmo.angDist(zss[0]) << std::endl;
  
  // This finds the 3D bounding box of the simulation particles
  Point_3d<float> Xmax,Xmin;
  halos[0].getBoundingBox(Xmin,Xmax);
  // gets the bounding box of the simulation particles
  for(int i=1 ; i < halos.size() ; ++i){
    Point_3d<float> p1,p2;
    halos[i].getBoundingBox(p1,p2);
    if(Xmax[0] < p2[0]) Xmax[0] = p2[0];
    if(Xmax[1] < p2[1]) Xmax[1] = p2[1];
    if(Xmax[2] < p2[2]) Xmax[2] = p2[2];
    
    if(Xmin[0] > p1[0]) Xmin[0] = p1[0];
    if(Xmin[1] > p1[1]) Xmin[1] = p1[1];
    if(Xmin[2] > p1[2]) Xmin[2] = p1[2];
  }

  std::cout << "Xmin = " << Xmin << std::endl;
  std::cout << "Xmax = " << Xmax << std::endl;

  Point_3d<float> center3d = (Xmax + Xmin)/2;
 
  // position of the densest particle
  //Point_3d<double> center3d = halomaker.densest_particle();
  
    Lens lens(&seed,zss[0],cosmo);
  
    // transfer the halos to the lens
    
    for(auto &h : halos){
      h.rotate(theta[0],theta[1],target);  // rotate 
      lens.moveinMainHalo(h, true);  // add the halos to the lens
      // Warning: this is not a copy of the LensHalo, but a move 
      // so the LensHalo in the halos vector is now empty.
      // By moving the halos into the lens one avoids
      // copying what might be a large amount of data.
      // The LensHalo's in the lens can be accessed with 
      // lens.getMainHalo(i) or lens.getMainHalos<LensHaloParticles>(i).
    }
  
  float range = (Xmax[0]-Xmin[0]) / Dl / 15.;
  //Point_2d angular_center(center3d[0]/Dl+range/2
  //  ,center3d[1]/Dl+range/2);
  
  // The center of the simulation is at (0,0) in the angular coordinates.
  // This can be changed with lens.getMainHalos<LensHaloParticles>(i).setTheta()
  Point_2d angular_center(0,0);

  std::cout << "angular center = " << angular_center << std::endl;
  std::cout << "range = " << range/arcminTOradians << " arcmin" << std::endl;
  logparams ("angular range of GridMap (arcmin)",range/arcminTOradians);

  std::string filename_tmp;
  std::vector<ImageFinding::CriticalCurve> critcurves;
  for(double zs : zss){
    
    filename_tmp = output_dir + fileprefix + "Zs"  + to_string(zs) + "prj" + to_string(projection);
    
    lens.ResetSourcePlane(zs);
    
    // here we shoot the rays in a grid
    start_time = time(nullptr);
    std::cout << "Making grid ..." ;
    GridMap grid(&lens,Npix,angular_center.x,range);

    end_time = time(nullptr);
    std::cout << "time to make grid " << difftime(end_time,start_time)/60 << " minutes" << std::endl;
    logparams("time to make grid",difftime(end_time,start_time) );
    
    // output maps of the lensing quantities
    if(do_maps){
      grid.writeFits<float>(LensingVariable::KAPPA,filename_tmp+".kappa.fits");
      grid.writeFits<float>(LensingVariable::ALPHA1,filename_tmp+".alpha1.fits");
      grid.writeFits<float>(LensingVariable::ALPHA2,filename_tmp+".alpha2.fits");
      grid.writeFits<float>(LensingVariable::INVMAG,filename_tmp+".invmag.fits");
      grid.writeFits<float>(LensingVariable::GAMMA1,filename_tmp+".gamma1.fits");
      grid.writeFits<float>(LensingVariable::GAMMA2,filename_tmp+".gamma2.fits");
    }
    
    // now we will find the critical curves
    int Ncrits;
    std::cout << std::endl << "Finding critical curves ..." ;
    start_time = time(nullptr);
    ImageFinding::find_crit(lens,grid,critcurves);
    end_time = time(nullptr);
    logparams("time to find critical curves",difftime(end_time,start_time) );
    
    std::cout << "found " << critcurves.size() << " critical curves" << std::endl;
    
    // print critical/caustic curve information to a file
    ImageFinding::printCriticalCurves(filename_tmp+"INFO",critcurves);
    
  }
  
  time_t t2 = time(nullptr);
  
  std::cout << "time = " << std::difftime(t2,t)/60 << " min"
  << std::endl;
  logparams("time (min)",std::difftime(t2,t)/60);
}



