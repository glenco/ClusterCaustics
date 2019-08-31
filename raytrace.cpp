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
    exit(0);
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
    
    COSMOLOGY cosmo(Planck18);
    
    cosmo.setOmega_matter(0.24,true);
    cosmo.sethubble(0.72);
    
    const int Npix =  2*500;
    const int Nsmooth = 30;
    const bool los = false;  // line-of-sight structure
    const bool cluster_on = true;
    const bool do_maps = false;
    long seed = -11920;

    const int projection = 1; // 1,2 or 3
    
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
    
    //std::string filename = "DataFiles/snap_058_centered.txt";
    std::string filename = "DataFiles/snap_058";
    //std::string filename = "DataFiles/snap_058_100000.txt";
    MakeParticleLenses halomaker(
                                 filename
                                 ,gadget2,Nsmooth,false
                                 );
    
    //std::string filename = "DataFiles/snap_69";
    
    //MakeParticleLenses halomaker(
    //                             filename
    //                             ,gadget2,Nsmooth,false
    //                             );

//    const double zl = 0.5294;
    const double zl = halomaker.getZoriginal();

    DLSDS func(cosmo,zl);
    //double zs = Utilities::bisection_search<DLSDS,double>(func
    //                                                      ,0.5,zl,10 ,0.001);
    
    
    double zs = 3;
    double Dl = cosmo.angDist(zl);
    
    std::cout << "zs = " << zs << std::endl;
    std::cout << "Dls/Ds = " << cosmo.angDist(zl, zs)/cosmo.angDist(zs) << std::endl;
    
    Point_3d Xmax,Xmin;
    halomaker.getBoundingBox(Xmin, Xmax);
    std::cout << "Xmin = " << Xmin << std::endl;
    std::cout << "Xmax = " << Xmax << std::endl;
    //Point_3d center3d = (Xmax + Xmin)/2;
    Point_3d center3d = halomaker.densest_particle();
    
    Point_3d target(500689.656250 , 498223.187500, 494912.468750);
    target *= 1.0e-3/(1+zl);///cosmo.gethubble();
 
    halomaker.Recenter(target);
    target *= 0;
    
    std::cout << center3d << std::endl;
    Point_2d center(0,0);
    std::cout << center << std::endl;
    
    // cut out a cylinder, could also do a ball
    //halomaker.cylindricalCut(center,(Xmax[0]-Xmin[0])/2);
    halomaker.radialCut(target,10);
    
    //long seed = 88277394;
    Lens lens(&seed,zs);
    
    //range /= cosmo.gethubble();
    //double range = (Xmax[0]-Xmin[0])*1.05/cosmo.gethubble()/Dl; // angular range of simulation
    center *= 1.0/cosmo.gethubble()/Dl; // convert to angular coordinates
    
    // 400 arcsecond range
    //double range = 4*300*arcsecTOradians;
    double range = 400*arcsecTOradians;

    std::cout << "area on sky " << range*range/arcsecTOradians/arcsecTOradians
    << " arcsec^2" << std::endl;
    
    filename = filename + ".sph" + to_string(Npix) + "x" + to_string(Npix) +
    "S" + to_string(Nsmooth) + "Zl"  + to_string(zl) + "Zs"  + to_string(zs) + "prj" + to_string(projection);

    if(cluster_on){
      halomaker.CreateHalos(cosmo,zl);
      for(auto h : halomaker.halos){
        h->rotate(theta);
        lens.moveinMainHalo(*h, true);
      }
    }else{
      filename = filename + "NoCL";
    }
    
    
    if(los){
      std::cout << "Generate LOS structures." << std::endl;
      lens.GenerateFieldHalos(1.0e11, ShethTormen
                              ,PI*range*range/2/degreesTOradians/degreesTOradians
                              ,20,nfw_lens,nsie_gal,2);
      filename = filename + "LOS10";
    }
    
    std::cout << "Making grid ..." ;
    Grid grid(&lens,Npix,center.x,range);
    std::vector<ImageFinding::CriticalCurve> critcurves;
    int Ncrits;
    
    if(do_maps){
      grid.writeFits(1,KAPPA,"!" + filename);
      grid.writeFits(1,ALPHA1,"!" + filename);
      grid.writeFits(1,ALPHA2,"!" + filename);
      grid.writeFits(1,INVMAG,"!" + filename);
    }
    std::cout << std::endl << "Finding critical curves ..." ;
    ImageFinding::find_crit(&lens,&grid,critcurves,&Ncrits,range/Npix/2);
    
    std::cout << "found " << Ncrits << " critical curves" << std::endl;
    ImageFinding::printCriticalCurves(filename,critcurves);
    
    if(critcurves.size() > 0){
      PixelMap map = ImageFinding::mapCausticCurves(critcurves,512*2*2);
      map.printFITS("!" + filename + "caust.fits");

      map = ImageFinding::mapCriticalCurves(critcurves,512*2*2);
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
      
      // ***************************************************
      //     image perturbations
      // ***************************************************
      
      const int Nsources = 500;
      const double r_source = 0.05; // in arcsec
      
      /// find largest critical curve
      area = critcurves[0].critical_area;
      ImageFinding::CriticalCurve &crit = critcurves[0];
      for(auto &a : critcurves){
        if(a.critical_area > area){
          area = a.critical_area;
          crit = a;
        }
      }
      
      vector<Point_2d> ys;
      Utilities::RandomNumbers_NR ran(seed);
      crit.RandomSourceWithinCaustic(Nsources, ys, ran);

      std::ofstream file_def(filename + "def.csv");

      std::vector<std::vector<Point_2d> > image_pos(Nsources);
      for(int j = 0; j < ys.size() ; ++j ){
        int Nimages;
        size_t Nimagepoints;
        std::vector<ImageInfo> imageinfo;
        ImageFinding::find_images_kist(&lens,ys[j].x, r_source * arcsecTOradians, &grid ,&Nimages,imageinfo,&Nimagepoints,
                                    100 * arcsecTOradians,true,1);
        
        for(int i=0; i < Nimages ; ++i){
          Point_2d p(imageinfo[i].centroid);
          cout << j << "  " << i << "  " << p << "  " << imageinfo[i].area/r_source/r_source
          /arcsecTOradians/arcsecTOradians/PI << endl;
          image_pos[j].push_back(p);
        }
      }
      
      lens.GenerateFieldHalos(1.0e11, ShethTormen
                              ,PI*range*range/2/degreesTOradians/degreesTOradians
                              ,20,nfw_lens,nsie_gal,2);
      Grid grid2(&lens,Npix,center.x,range);
      for(int j = 0; j < ys.size() ; ++j ){
        int Nimages;
        size_t Nimagepoints;
        std::vector<ImageInfo> imageinfo;
        ImageFinding::find_images_kist(&lens,ys[j].x, r_source * arcsecTOradians, &grid2 ,&Nimages
                                       ,imageinfo,&Nimagepoints, 100 * arcsecTOradians,true,1);
        
        for(auto im : imageinfo) cout << "     " << im.centroid[0] << " " << im.centroid[1] << " - " <<
          im.area/r_source/r_source/arcsecTOradians/arcsecTOradians/PI << endl;
        
        if(Nimages > image_pos[j].size() ){
          
          int k = 0;
          for(Point_2d p : image_pos[j] ){
            Point_2d delta(1.0e6,1.0e6);
            double length2 = delta.length_sqr();
            int imax = 0;
            for(int i=0; i < Nimages ; ++i){
              Point_2d p2(imageinfo[i].centroid);
              if( length2 > (p-p2).length_sqr() ){
                delta = p - p2;
                length2= delta.length_sqr();
                imax = i;
              }
            }
            cout << j << "  " << k << "  " << p << " " << delta << "  " << 1.0/imageinfo[imax].aveInvMag() << endl;
            file_def << j << "  " << k << "  " << p << " " << delta << "  " << 1.0/imageinfo[imax].aveInvMag() << endl;
            ++k;
          }
          
        }else{
          
          for(int i=0; i < Nimages ; ++i){
            Point_2d p2(imageinfo[i].centroid);
            Point_2d delta(1.0e6,1.0e6);
            double length2 = delta.length_sqr();
            int k = 0, kk = 0;
            for(Point_2d p : image_pos[j] ){
              if( length2 > (p-p2).length_sqr() ){
                delta = p-p2;
                length2 = delta.length_sqr();
                k = kk;
              }
              ++kk;
            }
            cout << j << "  " << i << " " << image_pos[j][k] << " " << delta << " " << 1.0/imageinfo[i].aveInvMag() << endl;
            file_def << j << " " << i << " " << image_pos[j][k] << " " << delta << " " << 1.0/imageinfo[i].aveInvMag()
            << endl;
          }
          
        }
       }
    
      file_def.close();
    }
    
    
    time_t t2 = time(nullptr);
    
    std::cout << "time = " << std::difftime(t2,t)/60 << " min"
    << std::endl;
  }

}



