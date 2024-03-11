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



//inline bool exists (const std::string& name) {
//  struct stat buffer;
//  return (stat (name.c_str(), &buffer) == 0);
//}

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
  
  
  time_t t;
  time(&t);
  
  COSMOLOGY cosmo(CosmoParamSet::Planck18);
  
  cosmo.setOmega_matter(0.24,true);
  cosmo.sethubble(0.72);
  
  const int Npix =  1000;// number of pixels in initial grid of rays
  const int Nsmooth = 30;// nearest neighbors smoothing scale
  const bool los_on = false;  // line-of-sight structure
  
  const bool cluster_on = true;  // toggle including cluster simulation
  const bool do_maps = true;    // map output maps or no
  //const std::vector<double> zss = {0.6,0.7,0.8,0.9,1.0,1.5,2.0};         // source redshift
  //const std::vector<double> zss = {2.0,1.5,1.0,0.9,0.8,0.7,0.6,0.5};         // source redshift
  const std::vector<double> zss = {2.0};         // source redshift
   //const std::vector<double> zss = {1.0};         // source redshift
  const bool no_images_perturb = true;
  
  const double range = 20*arcminTOradians;

  //long seed = -11920;
  long seed = time(&t);
  
  //which projection of simulation
  const int projection = 3; // 1,2 or 3
  
  //std::string filename = "DataFiles/snap_058_centered.txt";
  const std::string inputfilename = "DataFiles/snap_058";
  const std::string output_dir = "NewOutputMaps/";
  const std::string fileprefix = "snap_058";
  
  {
    std::ofstream logfile(output_dir + "logfile.txt");
    
    logfile << "Nsmooth : " << Nsmooth << endl;
    logfile << "projection : " << projection << endl;
    logfile << "zs : ";
    for(double zs : zss) logfile << zs << " ";
    logfile << endl;
    logfile << "los_on : " << los_on << endl;
    logfile << "cluster_on : " << cluster_on << endl;
    logfile << "do_maps : " << do_maps << endl;
    
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
  
  // this object is used to interface with the simulation file
  // , do the smoothing and make LensHalos for each type of
  // particle
  
  MakeParticleLenses halomaker(
                               inputfilename
                               ,SimFileFormat::gadget2   // the sim file format
                               ,Nsmooth
                               ,false     // takes particle type into account
                               );
  
  // gets the original redshift that was in the sim file
  const double zl = halomaker.getZoriginal();
  
  //DLSDS func(cosmo,zl);
  //double zs = Utilities::bisection_search<DLSDS,double>(func
  //                                                      ,0.5,zl,10 ,0.001);
  
  // angular size distance to the lens redshift
  double Dl = cosmo.angDist(zl);
  
  std::cout << "zs = " << zss[0] << std::endl;
  std::cout << "zl = " << zl << std::endl;
  std::cout << "Dls/Ds = " << cosmo.angDist(zl,zss[0])/cosmo.angDist(zss[0]) << std::endl;
  
  Point_3d<double> Xmax,Xmin;
  // gets the bounding box of the simulation particles
  halomaker.getBoundingBox(Xmin, Xmax);
  std::cout << "Xmin = " << Xmin << std::endl;
  std::cout << "Xmax = " << Xmax << std::endl;
  
  //Point_3d center3d = (Xmax + Xmin)/2;
  // position of the densest particle
  Point_3d<double> center3d = halomaker.densest_particle();
  
  // set position around which you want to extract particles
  // i.e. the center of the lens
  Point_3d<double> target(500689.656250 , 498223.187500, 494912.468750);
  target *= 1.0e-3/(1+zl);///cosmo.gethubble();
  
  // recenters simulation to target
  halomaker.Recenter(target);
  target *= 0;
  
  std::cout << center3d << std::endl;
  Point_2d center(0,0);
  std::cout << center << std::endl;
  
  // cut out a cylinder
  //halomaker.cylindricalCut(center,(Xmax[0]-Xmin[0])/2);
  
  //cut out a ball centered on target
  halomaker.radialCut(target,20);
  
  //range /= cosmo.gethubble();
  //double range = (Xmax[0]-Xmin[0])*1.05/cosmo.gethubble()/Dl; // angular range of simulation
  center *= 1.0/cosmo.gethubble()/Dl; // convert to angular coordinates
  
  // 400 arcsecond range
  //double range = 4*300*arcsecTOradians;
  
  std::cout << "area on sky " << range*range/arcsecTOradians/arcsecTOradians
  << " arcsec^2" << std::endl;
  
  Lens lens(&seed,zss[0],cosmo);
  
  if(cluster_on){
    // create the LensHalos, there will be one LensHalo
    // for each type of particle in this case
    halomaker.CreateHalos(cosmo,zl);
    // transfer the halos to the lens
    
    for(auto h : halomaker.halos){
      h->rotate(theta);
      lens.moveinMainHalo(*h, true);
      // By moving the halos into the lens one avoids
      // copying what might be a large amount of data.
      // The LensHalo in halomaker is now void.
    }
  }
  
  // this generates random halos along the line-of-sight
  if(los_on){
    std::cout << "Generate LOS structures." << std::endl;
    lens.GenerateFieldHalos(1.0e11,MassFuncType::ShethTormen
                            ,PI*range*range/2/degreesTOradians/degreesTOradians
                            ,20,LensHaloType::nfw_lens,GalaxyLensHaloType::nsie_gal,2);
  }
  
  std::string filename_tmp;
  std::vector<ImageFinding::CriticalCurve> critcurves;
  for(double zs : zss){
    
    filename_tmp = output_dir + fileprefix + ".sph" + to_string(Npix) + "x" + to_string(Npix) +
    "S" + to_string(Nsmooth) + "Zl"  + to_string(zl) + "Zs"  + to_string(zs) + "prj" + to_string(projection);
    
    if(!cluster_on) filename_tmp = filename_tmp + "NoCL";
    if(los_on) filename_tmp = filename_tmp + "LOS10";
    
    lens.ResetSourcePlane(zs);
    
    // here we shoot the rays in an initial grid
    
    std::cout << "Making grid ..." ;
    Grid grid(&lens,Npix,center.x,range);
    
    // output maps of the lensing quantities
    if(do_maps){
      grid.writeFits(1,LensingVariable::KAPPA,"!" + filename_tmp);
      //grid.writeFits(1,LensingVariable::ALPHA1,"!" + filename_tmp);
      //grid.writeFits(1,LensingVariable::ALPHA2,"!" + filename_tmp);
      grid.writeFits(1,LensingVariable::INVMAG,"!" + filename_tmp);
      grid.writeFits(1,LensingVariable::GAMMA1,"!" + filename_tmp);
      grid.writeFits(1,LensingVariable::GAMMA2,"!" + filename_tmp);
      PixelMap magmap = grid.writePixelMap(LensingVariable::INVMAG);
      
      for(size_t i=0 ; i<magmap.size() ; ++i) magmap[i] = 1/magmap[i];
      magmap.printFITS("!" + filename_tmp + ".mag.fits");
      
      PixelMap kappamap = grid.writePixelMap(LensingVariable::KAPPA);
      PixelMap gammamap = grid.writePixelMap(LensingVariable::GAMMA1);
      
      for(size_t i=0 ; i<kappamap.size() ; ++i) gammamap[i] /= (1+kappamap[i]);
      gammamap.printFITS("!" + filename_tmp + ".g1.fits");
      
      gammamap = grid.writePixelMap(LensingVariable::GAMMA2);
      
      for(size_t i=0 ; i<kappamap.size() ; ++i) gammamap[i] /= (1+kappamap[i]);
      gammamap.printFITS("!" + filename_tmp + ".g2.fits");

    }
    
    // now we will find the critical curves
    int Ncrits;
    std::cout << std::endl << "Finding critical curves ..." ;
    ImageFinding::find_crit(&lens,&grid,critcurves,&Ncrits,range/Npix/2);
    
    std::cout << "found " << Ncrits << " critical curves" << std::endl;
    
    // print critical/caustic curve information to a file
    ImageFinding::printCriticalCurves(filename_tmp+"INFO",critcurves);
    
    // make maps of the critical curves
    if(critcurves.size() > 0){
      //PixelMap map = ImageFinding::mapCausticCurves(critcurves,512*2*2);
      //map.printFITS("!" + filename_tmp + "caust.fits");
      
      //map = ImageFinding::mapCriticalCurves(critcurves,512*2*2);
      //map.printFITS("!" + filename_tmp + "crit.fits");
      
      // output points along caustic curve
      std::ofstream file(filename_tmp + "caustic.csv");
      double area = 0;
      int j=0;
      for(int i = 0 ; i < critcurves.size() ; ++i){
        if(critcurves[i].critical_area > area){
          area = critcurves[i].critical_area;
          for(RAY p : critcurves[i].critcurve){
            file << j << "," << p << std::endl;
          }
          ++j;
        }
      }
      file.close();
      
    }
  }

  if(no_images_perturb) exit(1);

  // ***************************************************
  //     image perturbations
  //
  //   Here I am generating random LOS objects and seeing
  //   how the image positions change.
  // ***************************************************
  
  const int Nsources = 500;
  const double r_source = 0.05; // in arcsec
  const bool fixed_image = true; // fixes one image to the same place when the the LSS is present or not
  const string deflection_sufix= "fixdef.csv";
  
  /// find largest critical curve
  double area = critcurves[0].critical_area;
  ImageFinding::CriticalCurve &crit = critcurves[0];
  for(auto &a : critcurves){
    if(a.critical_area > area){
      area = a.critical_area;
      crit = a;
    }
  }
  
  vector<Point_2d> ys;
  vector<Point_2d> index_image;
  
  Utilities::RandomNumbers_NR ran(seed);
  crit.RandomSourceWithinCaustic(Nsources, ys, ran);
  
  Grid grid(&lens,Npix,center.x,range);
  {
    std::ofstream file_def(filename_tmp + deflection_sufix);
    file_def << "lens image x_image y_image delta_x delta_y mag mag_sign same_number" << endl;
    
    std::vector<std::vector<Point_2d> > image_pos(Nsources);
    for(int j = 0; j < ys.size() ; ++j ){
      int Nimages;
      size_t Nimagepoints;
      std::vector<ImageInfo> imageinfo;
      ImageFinding::find_images_kist(&lens,ys[j].x, r_source * arcsecTOradians, &grid ,&Nimages,imageinfo,&Nimagepoints,
                                     100 * arcsecTOradians,true,1);
      
      index_image.push_back(imageinfo[0].centroid);
      for(int i=0; i < Nimages ; ++i){
        Point_2d p(imageinfo[i].centroid);
        cout << j << "  " << i << "  " << p << "  " << imageinfo[i].area/r_source/r_source
        /arcsecTOradians/arcsecTOradians/PI << endl;
        image_pos[j].push_back(p);
      }
    }
    
    lens.GenerateFieldHalos(1.0e11,MassFuncType::ShethTormen
                            ,PI*range*range/2/degreesTOradians/degreesTOradians
                            ,20,LensHaloType::nfw_lens,GalaxyLensHaloType::nsie_gal,2);
    
    LinkedPoint point;
    if(fixed_image){
      int i=0;
      for(Point_2d p : index_image){
        point.x[0] = p[0]; point.x[1] = p[1];
        lens.rayshooterInternal(1,&point);
        ys[i][0] = point.image->x[0];
        ys[i++][1] = point.image->x[1];
      }
    }
    
    
    Grid grid2(&lens,Npix,center.x,range);
    for(int j = 0; j < ys.size() ; ++j ){
      int Nimages;
      size_t Nimagepoints;
      std::vector<ImageInfo> imageinfo;
      ImageFinding::find_images_kist(&lens,ys[j].x, r_source * arcsecTOradians, &grid2 ,&Nimages
                                     ,imageinfo,&Nimagepoints, 100 * arcsecTOradians,true,1);
      
      for(auto im : imageinfo) cout << "     " << im.centroid[0] << " " << im.centroid[1] << " - " <<
        im.area/r_source/r_source/arcsecTOradians/arcsecTOradians/PI << endl;
      
      bool same_number = Nimages == image_pos[j].size();
      
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
          cout << j << "  " << k << "  " << p << " " << delta << " "
          << imageinfo[imax].area/r_source/r_source/arcsecTOradians/arcsecTOradians/PI << " "
          << (imageinfo[imax].aveInvMag() > 0)*2 - 1 << " "
          << same_number
          << endl;
          
          file_def << j << " " << k << " " << p << " " << delta << " "
          << imageinfo[imax].area/r_source/r_source/arcsecTOradians/arcsecTOradians/PI << " "
          << (imageinfo[imax].aveInvMag() > 0)*2 - 1 << " "
          << same_number
          << endl;
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
          cout << j << " " << i << " " << image_pos[j][k] << " " << delta << " "
          << imageinfo[i].area /r_source/r_source/arcsecTOradians/arcsecTOradians/PI << " "
          << (imageinfo[i].aveInvMag() > 0)*2 - 1 << " "
          << same_number
          << endl;
          
          file_def << j << " " << i << " " << image_pos[j][k] << " " << delta << " "
          << imageinfo[i].area /r_source/r_source/arcsecTOradians/arcsecTOradians/PI << " "
          << (imageinfo[i].aveInvMag() > 0)*2 - 1 << " "
          << same_number
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



