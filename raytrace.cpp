/* ParticleExample
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
#include "particle_halo.h"
#include "utilities.h"
#include "concave_hull.h"
//#include "gadget.hh"

using namespace std;



inline bool exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

struct DLSDS{
  DLSDS(COSMOLOGY &c,double zl):zl(zl),cosmo(c){

   // Dl = cosmo.angDist(zl);
  };
  
  double zl;
  COSMOLOGY &cosmo;
  
  double operator()(double zs){return cosmo.angDist(zl,zs)/cosmo.angDist(zs);}
};

int main(int arg,char **argv){
  
  /*{  // Test of critical curve ordering
    std::vector<double> data;
    int ncol,nrows;
    std::string filename = "DataFiles/snap_058_centered.txt.cy2049x2049S60Zl0.506000caustic.csv";
  
    Utilities::IO::ReadASCII(data, filename, ncol, nrows);
  
    std::cout << ncol << " " << nrows << std::endl;
    
    std::vector<Point_2d> points(nrows);
    for(int i = 0 ; i < nrows ; ++i){
      points[i][0] = data[ ncol*i];
      points[i][1] = data[ ncol*i + 1];
    }
    
    double scale = 10*sqrt( pow(points[0][0] - points[1][0],2) + pow(points[0][1] - points[1][1],2) );
    std::vector<Point_2d> hull;
    Utilities::concave(points,hull,scale);
    //Utilities::convex_hull(points,hull);
    //std::vector<Point_2d> hull = Utilities::concave2(points, scale);
    
    std::ofstream file("caustic_test.csv");
    file << "x,y" << std::endl;
    for(auto p : hull){
      file << p[0] << "," << p[1] << std::endl;
    }
    file.close();
    
  }*/
  
  /*{
    std::vector<ParticleType<> > data;
    
    GadgetFile<ParticleType<float> > gadget_file("DataFiles/snap_069",data);
    gadget_file.openFile();
    gadget_file.readBlock("POS");
    gadget_file.readBlock("MASS");

    std::cout << data[10][0] << " " << data[10][1] << " " << data[10][2] << std::endl;
    std::cout << data[10].mass() << " " << data[10].size() << std::endl;
    exit(0);
  }*/
  
  
  {
   
    time_t t;
    time(&t);
    
    COSMOLOGY cosmo(Planck);
    
    const double zl = 0.506;
    const int Npix =  2*2049;
    const int Nsmooth = 60;
    const bool los =false;
    
    DLSDS func(cosmo,zl);
    double zs = Utilities::bisection_search<DLSDS,double>(func
                                ,0.5,zl,10,0.001);

    double Dl = cosmo.coorDist(zl);

    std::cout << "zs = " << zs << std::endl;
    std::cout << "Dls/Ds = " << cosmo.angDist(zl, zs)/cosmo.angDist(zs) << std::endl;
    
    std::string filename = "DataFiles/snap_058_centered.txt";
    //std::string filename = "DataFiles/snap_058_short.txt";
    MakeParticleLenses halomaker(
                                 filename
                                 ,csv4,Nsmooth,false
                                 );
    
    //std::string filename = "DataFiles/snap_69";
    
    //MakeParticleLenses halomaker(
    //                             filename
    //                             ,gadget2,Nsmooth,false
    //                             );

    
    Point_3d Xmax,Xmin;
    halomaker.getBoundingBox(Xmin, Xmax);
    
    Point_3d center3d = (Xmax + Xmin)/2;
    Point_2d center;
    center[0] = center3d[0];
    center[1] = center3d[1];
    // cut out a cylinder, could also do a ball
    halomaker.cylindricalCut(center,(Xmax[0]-Xmin[0])/2);
    
    //long seed = 88277394;
    long seed = -11920;
    Lens lens(&seed,zs);
    
    double range = (Xmax[0]-Xmin[0])*1.05*cosmo.gethubble()/Dl; // angular range of simulation
    center *= cosmo.gethubble()/Dl; // convert to angular coordinates
    
    std::cout << "area on sky " << range*range/arcsecTOradians/arcsecTOradians
    << " arcsec^2" << std::endl;
    
    halomaker.CreateHalos(cosmo,zl);
    for(auto h : halomaker.halos){
      lens.insertMainHalo(h,zl, true);
    }
    
    filename = filename + ".cy" + to_string(Npix) + "x" + to_string(Npix) +
    "S" + to_string(Nsmooth) + "Zl"  + to_string(zl);
    
    if(los){
      lens.GenerateFieldHalos(1.0e11, ShethTormen
                              ,PI*range*range/4/degreesTOradians/degreesTOradians
                              ,20,nfw_lens,nsie_gal,2);
      filename = filename + "LOSg";
    }
    
    Grid grid(&lens,Npix,center.x,range);
    std::vector<ImageFinding::CriticalCurve> critcurves;
    int Ncrits;
    
    ImageFinding::find_crit(&lens,&grid,critcurves,&Ncrits,range/Npix);
    
    std::cout << "found " << Ncrits << " critical curves" << std::endl;
    ImageFinding::printCriticalCurves(filename,critcurves);
    
    grid.writeFits(1,KAPPA,"!" + filename);
    grid.writeFits(1,ALPHA1,"!" + filename);
    grid.writeFits(1,ALPHA2,"!" + filename);
    grid.writeFits(1,INVMAG,"!" + filename);
    

    PixelMap map = ImageFinding::mapCausticCurves(critcurves,512*2);
    map.printFITS("!" + filename + "caust.fits");

    map = ImageFinding::mapCriticalCurves(critcurves,512*2);
    map.printFITS("!" + filename + "crit.fits");

    double area = 0;
    int j=0;
    for(int i = 0 ; i < critcurves.size() ; ++i){
      if(critcurves[i].critical_area > area){
        area = critcurves[i].critical_area;
        j = i;
      }
    }
    
    std::ofstream file(filename + "caustic.csv");
    for(Point_2d p : critcurves[j].critical_curve){
      file << p << std::endl;
    }
    file.close();
 
    time_t t2 = time(nullptr);
    
    std::cout << "time = " << std::difftime(t2,t)/60 << " min"
    << std::endl;
  }

}



