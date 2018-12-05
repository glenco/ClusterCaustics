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
    const int Npix =  2049;
    const int Nsmooth = 60;
    const bool los = true;
    
    DLSDS func(cosmo,zl);
    double zs = Utilities::bisection_search<DLSDS,double>(func
                                ,0.5,0,5,0.001);

    double Dl = cosmo.coorDist(zl);

    std::cout << "zs = " << zs << std::endl;
    
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
      lens.GenerateFieldHalos(1.0e11, ShethTormen,2*PI*range*range/4/degreesTOradians/degreesTOradians
                              ,20,nfw_lens,nsie_gal);
      filename = filename + "LOSg";
    }
    
    Grid grid(&lens,Npix,center.x,range);
    std::vector<ImageFinding::CriticalCurve> critcurves;
    int Ncrits;
    
    ImageFinding::find_crit(&lens,&grid,critcurves,&Ncrits,range/Npix/2);
    
    std::cout << "found " << Ncrits << " critical curves" << std::endl;
    ImageFinding::printCriticalCurves(filename,critcurves);
    
    grid.writeFits(1,KAPPA,"!" + filename);
    grid.writeFits(1,ALPHA1,"!" + filename);
    grid.writeFits(1,ALPHA2,"!" + filename);
    grid.writeFits(1,INVMAG,"!" + filename);
    
    /*
    GridMap gridmap(&lens,Npix,center.x,range);
    
    PixelMap pmap = gridmap.writePixelMapUniform(KAPPA);
    pmap.printFITS("!" + filename + ".kappa.fits");

    pmap = gridmap.writePixelMapUniform(ALPHA1);
    pmap *= 1.0/arcsecTOradians;
    pmap.printFITS("!" + filename + ".alpha1.fits");

    pmap = gridmap.writePixelMapUniform(ALPHA2);
    pmap *= 1.0/arcsecTOradians;
    pmap.printFITS("!" + filename + ".alpha2.fits");
    
    pmap = gridmap.writePixelMapUniform(INVMAG);
    pmap.printFITS("!" + filename + ".invmag.fits");
    */

    time_t t2 = time(nullptr);
    
    std::cout << "time = " << std::difftime(t2,t)/60 << " min"
    << std::endl;
    exit(0);
  }
  
  /************************************************
   {
   
   PointList plist;
   PointList::iterator it;
   Point p;
   plist.InsertPointAfterCurrent(it,&p);
   for( it = plist.begin() ; it != plist.end() ; ++it){
   cout << " Hi " << endl;
   }
   
   PixelMap map_source("galaxy.fits",4.8e-7);
   //map_source.printFITS("source_copy.fits");
   
   
   std::cout << " read source from fits "  << std::endl;
   SourcePixelled source(map_source,3,1.0);
   source.setX(0.0,0.0);
   }
   //*************************************************/
  
  COSMOLOGY cosmo(Planck1yr);
  Point_2d rotation_vector;
  PosType zl=0.5;
  PosType zs = 2; //** redshift of source
  rotation_vector *= 0;
  if(arg < 2){
    std::cerr << " *** Needs command line arguments ***" << std::endl;
    std::cerr << "    1) ID number for" << std::endl;
    std::cerr << "    2) projection, axis that is the line of sight" << std::endl;
    std::cerr << "    3) smoothing number" << std::endl;
  }
  int Nsmooth=std::stoi(argv[3]);
  int projection=std::stoi(argv[2]);
  int subID = std::stoi(argv[1]);
  long seed = -28976391;
  
  /****************
  {
    const PixelMap map_source("galaxy.fits",4.8e-7);
    map_source.printFITS("!source_copy.fits");
    
    std::cout << " read source from fits "  << std::endl;
    SourcePixelled source(map_source,3,1);
    
    LensHaloDummy dhalo;
    Lens lens(&seed,5);
    lens.insertMainHalo(&dhalo,zl,true);
    GridMap gridmap(&lens,500,source.getCentroid(),map_source.getRangeX());
    gridmap.RefreshSurfaceBrightnesses(&source);
    PixelMap map_image = gridmap.getPixelMap(1.0);
    map_image.printFITS("!test_run_nolens.fits");
  }*/

  std::string filename(std::to_string(subID)+"_p"+std::to_string(projection)+"_ft");
  
  LensHaloParticles<ParticleType<float> > dm_halo("DataFiles/particles_main_"+std::to_string(subID)+"_dm.txt",ascii, zl, Nsmooth, cosmo, rotation_vector, false, false,0);
  LensHaloParticles<ParticleType<float> > gas_halo("DataFiles/particles_main_"+std::to_string(subID)+"_gas.txt",ascii, zl, Nsmooth, cosmo, rotation_vector, false, false,0);

  LensHaloParticles<ParticleType<float> > st_halo("DataFiles/particles_main_"+std::to_string(subID)+"_stars.txt",ascii, zl, Nsmooth, cosmo, rotation_vector, false, false,0);

  // Main halo w/o substructure^M
  //LensHaloParticles dm_halo("DataFiles/particles_main_"+std::to_string(subID)+"_dm.txt",zl,Nsmooth,cosmo,rotation_vector, false, true);
  //LensHaloParticles gas_halo("DataFiles/particles_main_"+std::to_string(subID)+"_gas.txt",zl,Nsmooth,cosmo,rotation_vector, false, true);
  //LensHaloParticles st_halo("DataFiles/particles_main_"+std::to_string(subID)+"_stars.txt",zl,Nsmooth,cosmo,rotation_vector, false, true);
  
  //LensHaloRealNSIE nsie_halo(1.0e12, 0.5, 250, 0, 0.5, 0, 0);

  
  // Main halo w/ substructure
  InputParams params("param_example");
  Lens lens(params,&seed);
  
  
  //double start,potential,src_st,src_ed;
  //start = clock();
  
  
  {
    // This is makes a SIE made out of particles for testing purposes.
    //  The MOCK data is saved to a file.
    //Utilities::RandomNumbers_NR ran(12083);
    //LensHaloParticles::makeSIE("test.dat",0.5,1.0e9,1.0e11,100,0.2,ran);
    //exit(0);
  }


  // Create Main halo w/o substructure
  
  std::cout << "Importing particles and calculating smoothing ..." << std::endl ;
  std::cout << "    This can take a while if it hasn't already been done on this data...." << std::endl ;
  
  // replaces object put in from parameter file
  
  Point_2d center,proj_center,theta;
  PosType Dl = cosmo.angDist(zl);
  
  // determine projection
  std::cout<<"subfindID  "<<subID<<"  projection  "<<projection<<endl;
  
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
 
    Point_3d cm_3d = dm_halo.CenterOfMass();
    std::cout << "center of mass before rotation" << std::endl;
    std::cout<<cm_3d[0]<<"  "<<cm_3d[1]<<"  "<<cm_3d[2]<<endl;

    // rotate the simulation
    dm_halo.rotate(theta);
    gas_halo.rotate(theta);
    st_halo.rotate(theta);
    
    cm_3d = dm_halo.CenterOfMass();
    std::cout << "center of mass after rotation" << std::endl;
    std::cout<<cm_3d[0]<<"  "<<cm_3d[1]<<"  "<<cm_3d[2]<<endl;

    proj_center[0] = cm_3d[0]/Dl;
    proj_center[1] = cm_3d[1]/Dl;
  }
  else{
    std::cout<<"projection input should be within 1-3!"<<endl;
    cout << "Exiting" << endl;
    exit(1);
  }

  //std::cout<<cm_3d[0]<<"  "<<cm_3d[1]<<"  "<<cm_3d[2]<<endl;
  
  // ******** Insert halos into lens ******

  //nsie_halo.setTheta(0,0);    // set halo position
  //lens.replaceMainHalo(&nsie_halo,0.4,true);
  
  lens.replaceMainHalo(&dm_halo,zl,true);
  lens.insertMainHalo(&gas_halo,zl,true);
  lens.insertMainHalo(&st_halo,zl,true);
  
  center = proj_center;
  
  // reset the redshift of the source plane from one already set in parameter file
  lens.ResetSourcePlane(zs,false);
  
  double range = 0.25*degreesTOradians/10; // range of grids in radians
 
  // create a Grid of rays.  You could also use GridMap if no further refinement
  // was required (see below).
  //center = {0,0};
  Grid grid(&lens,512,center.x,range/2);
  
  // output some maps

  grid.writeFits(1.0,KAPPA,"!particles" + filename + "_"+std::to_string(Nsmooth));
  grid.writeFits(1.0,INVMAG,"!particles"  + filename +"_"+ std::to_string(Nsmooth));
  grid.writeFits(1.0,DT,"!particles"  + filename +"_"+ std::to_string(Nsmooth));

  
  // find the critical curves
  std::vector<ImageFinding::CriticalCurve> crit_curve;
  int Ncriticals;
  std::cout << "Finding caustics..." << std::endl;
  ImageFinding::find_crit(&lens,&grid,crit_curve,&Ncriticals,0.01*arcsecTOradians);
  
  std::cout << "Number of caustics : "<< crit_curve.size() << std::endl;
  //std::cout << "Number of caustics : "<< Ncriticals << std::endl;

  if(crit_curve.size() == 0){
    cout << "Exiting" << endl;
    exit(1);
  }
  
  //*** plot caustic curves
  
  Point_2d caus_plot_center,crit_plot_center;
  double caus_range,crit_range;
  
  //*** print information about the critical curves that were found
  PosType rmax,rmin,rave;
  std::vector <PosType> tan_rmin,caus_rmax,crit_rmax;
  std::vector <int> tan_idx; // find main tangential caustic
  int idx = -1; // idx for main tangential caustic
  double tancaus_rmin;
  double caus_cx,caus_cy;
  
  if(crit_curve.size() > 0){
    std::string type;
    for(int i=0;i<crit_curve.size();++i){
      type = to_string( crit_curve[i].type );
      std::cout << "  " << i << " type " << to_string(crit_curve[i].type) << std::endl;
      crit_curve[i].CausticRadius(rmax,rmin,rave);
      
      if(to_string(crit_curve[i].type)=="tangential"){
        tan_rmin.push_back(rmin);
        tan_idx.push_back(i);
        std::cout << "Record  " << i << " type " << to_string(crit_curve[i].type) << std::endl;
      }
      caus_rmax.push_back(rmax);
      std::cout << "      caustic " << crit_curve[i].caustic_center << " | " << crit_curve[i].caustic_area << " " << rmax << " " << rmin << " " << rave << std::endl;
      // ------
      crit_curve[i].CriticalRadius(rmax,rmin,rave);
      std::cout << "      critical " << crit_curve[i].critical_center << " | " << crit_curve[i].critical_area << " " << rmax << " " << rmin << " " << rave << std::endl;
      crit_rmax.push_back(rmax);
    }
    
    // find main tangential caustics / largest critical index
    
    if(tan_idx.size()==0){
      std::cout << "It is not a standard lens!!!" << std::endl;
      std::cout << "No tangential caustic found." << std::endl;
      exit(0);
    }
    
    if(tan_idx.size()==1) idx=tan_idx[0];
    
    if(tan_idx.size()>1) {
      PosType max_rmin= *std::max_element(tan_rmin.begin(),tan_rmin.end());
      for(int i=0;i<tan_rmin.size();++i){
        if(tan_rmin[i]==max_rmin) idx=i;
      }
    }
    
    crit_curve[idx].CausticRadius(rmax,rmin,rave);
    tancaus_rmin=rmin;
    std::cout<<"The main tangential caustic is: "<<idx<<" rmin= "<< rmin<<endl;
    
    crit_curve[idx].CriticalRadius(rmax,rmin,rave);
    crit_range=1.1*2*rmax;
    crit_plot_center=crit_curve[idx].critical_center;
    std::cout<<"The largest critical is: "<<idx<<" rmax= "<< rmax<<endl;
    
    // get main caustic center
    caus_cx=crit_curve[idx].caustic_center[0];
    caus_cy=crit_curve[idx].caustic_center[1];
    
    // find largest caustic
    int causmax_id;
    //PosType max_rmax = *std::max_element(caus_rmax.begin(),caus_rmax.end());
    //for(int i=0;i<crit_curve.size();++i){
    //  if(caus_rmax[i]==max_rmax) causmax_id=i;
    //}
    PosType max_rmax = *std::max_element(crit_rmax.begin(),crit_rmax.end());
    for(int i=0;i<crit_curve.size();++i){
      if(crit_rmax[i]==max_rmax) causmax_id=i;
    }
    
    crit_curve[causmax_id].CausticRadius(rmax,rmin,rave);
    caus_range=1.1*2*rmax;
    caus_plot_center=crit_curve[causmax_id].caustic_center;
    std::cout<<"The largest caustic is: "<<causmax_id<<" rmax= "<< rmax<<endl;
  }
  
  //** write down main tangential caustic pts
  fstream caus_pt;
  caus_pt.open("caustic_"+filename+".txt",ios::out|ios::trunc);
  
  for(int i=0;i<crit_curve[idx].caustic_curve_outline.size();++i){
    caus_pt<< crit_curve[idx].caustic_curve_outline[i].x[0] <<" "<< crit_curve[idx].caustic_curve_outline[i].x[1] << endl;
  }
  caus_pt.close();
  
  Point_2d map_center;
  center=caus_plot_center;
  range=caus_range;
  PixelMap map(center.x,512*2,range/512/2);// caustics
  
  for(int i=0;i<crit_curve.size();++i){
    map.AddCurve(crit_curve[i].caustic_curve_intersecting, crit_curve[i].type+1);
  }
  map.printFITS("!caustics_"  + filename +"_"+std::to_string(Nsmooth) +"NN.fits");
  
  //*** plot critical curves
  center=crit_plot_center;
  range=crit_range;
  PixelMap map_crit(center.x,1024,range/1024);
  
  for(int i=0;i<crit_curve.size();++i) {
    map_crit.AddCurve(crit_curve[i].critical_curve, crit_curve[i].type);
  }
  
  
  map_crit.printFITS("!critical_"  + filename +"_" + std::to_string(Nsmooth) +"NN.fits");
  std::cout << " critical curve range = " << map_crit.getRangeX()/arcsecTOradians <<  " arc seconds " << map.getRangeX()/grid.getInitRange() << std::endl;
  
  PosType factor = 1;
  
  const PixelMap map_source("galaxy.fits",4.8e-7/100.);
  map_source.printFITS("!source_copy.fits");
  
  
  std::cout << " read source from fits "  << std::endl;
  SourcePixelled source(map_source,zs,factor);
  source.setTheta(crit_plot_center.x);
  
  
  {
    Point_2d xs,xim;
    source.getTheta(xs);
    map_crit.getCenter(xim);
    
    cout << " map range / source size = " << map_crit.getRangeX()/source.getSize() << endl;
    cout << " (source - map center)/map_crit.getRangeX()     =  " << (xs-xim).length()/map_crit.getRangeX() << endl;
    grid.RefreshSurfaceBrightnesses(&source);
    PixelMap map_image(map_crit);  // this way it will be the same size as the other maps
    map_image.Clean();
    //PixelMap map_image(center.x, 1024/2, 0.25*pi/180/10*2/1024);
    map_image.AddGridBrightness(grid);
    map_image.printFITS("!test_run_grid.fits");
  }
  //************************/
  
  GridMap gridmap(&lens,200,center.x,map_crit.getRangeX());
  gridmap.RefreshSurfaceBrightnesses(&source);
  
  
  //gridmap.getPixelMap("test_run.fits");
  PixelMap map_image = gridmap.getPixelMap(1.0);
  map_image.printFITS("!test_run.fits");
  
  //std::cout <<"gridmap done"<<std::endl;
  
  // A PixelMap is made for the image
  // the image is put on the pixel map
  /* PixelMap map_image = gridmap.getPixelMap(1.0);
   
   std::vector<ImageInfo> imageinfo;
   int Nimages;
   
   
   ImageFinding::map_images(&lens,&source,&grid,&Nimages,imageinfo,2.42e-5,2.42e-9,0,EachImage,false,false);
   map_image.printFITS("!image_"  + filename +"_"+std::to_string(Nsmooth) +".fits");
   */
  return 0;
}



