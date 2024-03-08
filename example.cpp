#include <slsimlib.h>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <thread>
#include <mutex>
#include "grid_maintenance.h"
#include "gridmap.h"

//std::ostream &operator<<(std::ostream &os, Point_2d const &p) {
// return os << p.x[0] << " " << p.x[1];
// }
 

int main(int arg,char **argv){

  
  //****************************************
  // **** random number seed, not actually used in this example, but
  // ***** needed for lens constructor
  //****************************************
  long seed = -1827675;
  double zeropoint = 45.0;  // arbitrary example

  //**** set the cosmology  ******
  COSMOLOGY cosmo(CosmoParamSet::Planck18);

  
  // construct some massive objects call "halos"
  // an NFW halo - mass,Rmax,redshift,concentration,axis ratio
  //               ,position angle,number of stars
  
  double z_lens = 0.3; // redshift of lens
  double z_source = 3; // redshift of source
  
  //LensHaloNFW halo_nfw(1.0e13,1.0,z_lens,4,1,0,0);
  //LensHaloNFW halo_nfw(1.0e13,1.0,z_lens,4,1,0,cosmo);
  //halo_nfw.setTheta(0,0); // set angular position
  // Truncated NonSingular Isothermal Ellipsoid - mass,redshift,sigma
  // , core radius,axis ratio,position angle,number of stars
  LensHaloTNSIE halo_nsie(1.0e13,z_lens, 300, 0, 0.5, 45, cosmo);
  halo_nsie.setTheta(0,0);

  // construct lens
  Lens lens(&seed,z_source,cosmo);
  
  // insert the LensHalos into the lens
  //lens.insertMainHalo(halo_nfw,true);
  lens.insertMainHalo(halo_nsie,true);
  
  std::cout << "Lens constructed" << std::endl;
  
  
  /*******************************
  
   Now that we have created a Lens we can shoot rays through it
   
  ********************************/
  //double center[] = {0.3*pi/180,-0.25*pi/180};
  double center[] = {0,0};  // center of grids
  //double range = 0.3*pi/180/20/2; // range of grids in radians
  double range = 15*arcsecTOradians; // range of grids in radians
  size_t Ninit = 1024/2/2; // the initial number of pixels to a side in the grid
  
  
  std::cout << " range = " << range/degreesTOradians << " degrees" << std::endl;
  /**
   Here a uniform grid is constructed that is not capable 
   of refinement for images, caustics, etc.  This is the 
   simplest way to make a map.
   In the construction the rays are shot through the lens 
   in parallel if the N_THREAD is set to multiple threads 
   in the GLAMER build.
  **/
  
    std::cout << "Constructing initial GridMap ..." << std::endl;
    GridMap gridmap(&lens,Ninit,center,range);
    std::cout << "constructed" << std::endl;

    /***
     make an images of kappa
     The center,range and pixel numbers do not have to match the grid, but 
     in this case they are set to match.
     ****/
    gridmap.writeFitsUniform(LensingVariable::KAPPA,"!kappa_map.fits");
    gridmap.writeFitsUniform(LensingVariable::INVMAG,"!invers_magnification_map.fits");


  /******************************************
   Now we are going to look for same caustics
   ******************************************/

  // The CriticalCurve class contains information about a caustic
  std::vector<ImageFinding::CriticalCurve> critcurves;
  
  std::cout << "Looking for critical curves ..." << std::endl;
  // Find all the critical curves
  ImageFinding::find_crit(lens,gridmap,critcurves);
  std::cout << critcurves.size() << " critical curves found." << std::endl;

  PosType Xrange[2]={0,0},Yrange[2]={0,0};
  PosType rmax,rmin,rave;
  for(auto crit : critcurves){
    crit.CriticalRadius(rmax,rmin,rave);
    std::cout << (int)(crit.type) << " " << crit.critical_area << " " << crit.caustic_area << std::endl;;
  }
  if(critcurves.size() == 0){
    std::cout << "No critical curves where found. Construct another lens." << std::endl;
  }else{
    Point_2d p1,p2;
    
    // find a box on the image plane that contains all of the critical curves
    for(auto crit : critcurves){
      crit.CritRange(p1,p2);
      
      Xrange[0] = MIN(Xrange[0],p1[0]);
      Xrange[1] = MAX(Xrange[1],p2[0]);
      
      Yrange[0] = MIN(Yrange[0],p1[1]);
      Yrange[1] = MAX(Yrange[1],p2[1]);
    }
    
    //***************************************
    // This part makes a fits image of the caustics and critical curves.
    // It is crude.  A better way to make a nice plot is to print the points
    // in the curve to a file and plot them with some other program.
    //***************************************
    
    // print critical curves and caustics to a csv file
    {
      std::ofstream caustics("caustics_and_crits.csv");
      int i=0;
      caustics << "# id - caustic number" << std::endl;
      caustics << "# type - type of caustic, radial,tangential,psudo" << std::endl;
      caustics << "# crit_x - x coordinate of point on critical curve (radians)" << std::endl;
      caustics << "# crit_y - y coordinate of point on critical curve (radians)" << std::endl;
      caustics << "# caust_x - x coordinate of point on caustic curve (radians)" << std::endl;
      caustics << "# caust_y - y coordinate of point on caustic curve (radians)" << std::endl;
      caustics << "id,type,crit_x,crit_y,caust_x,caust_y" << std::endl;
      for(ImageFinding::CriticalCurve &crit : critcurves){
        for(RAY &p : crit.critcurve){
          caustics << i << ",";
          caustics << crit.type << ",";
          caustics << p.x[0] << ",";
          caustics << p.x[1] << ",";
          caustics << p.y[0] << ",";
          caustics << p.y[1] << std::endl;
        }
        ++i;
      }
    }
    // *****************************************
    //  RAYs can also be shot individually or in groups in parallel through the lens.
    // They contain source in image position along with the local magnification matrix.
    // *****************************************
    
    //**************************
    // Now we are going to pick a caustic, put a source
    // in it and find its images
    //**************************
    
    std::vector<Point_2d> ys;  // source positions
    // this is a random number generator
    Utilities::RandomNumbers_NR rng(seed);
    
    // This finds random points within this caustic.
    // In this case it is just one point.
    int nc=0;  // pick the caustic with largest critical curve
    for (int i=1 ; i < critcurves.size() ; ++i) {
      if(critcurves[nc].critical_area > critcurves[i].critical_area) nc = i;
    }
    critcurves[nc].RandomSourceWithinCaustic(1,ys,rng); // 1 random point within caustic
    
    // ImageInfo is a class that contains information about images
    std::vector<ImageInfo> imageinfo;
    int Nimages;  // number of images that will be found
    size_t Nimagepoints;  // total number of points within the images
    
    std::vector<Point_2d> image_points;
    {
      std::vector<GridMap::Triangle> trs;
      gridmap.find_images(ys[0],image_points, trs);
    }
    std::cout << "Number of images: " << image_points.size() << std::endl;
    
    //******************************************************
    // This is a different way of finding images without going through a GridMap
    // It can find the images positions to whatever resolition is required.
    // In this case 1/5th the resolution of the GridMap.
    //******************************************************
    
    std::vector<RAY> rays = lens.find_images(ys[0], z_source, gridmap.getCenter(), gridmap.getXRange(), gridmap.getResolution()/5);
    
    for(RAY ray : rays){
      std::cout << ray << std::endl;
    }
    
    //*************************************************************
    //**** Let us put a more complicated source in the image
    //*************************************************************
    
    
    //*** find a source position within the tangential caustic
    Utilities::RandomNumbers_NR random(seed);   //*** random number generator
    std::vector<Point_2d> y;                    //*** vector for source positions
    critcurves[nc].RandomSourceWithinCaustic(2,y,random); // get random points within first caustic
    
    critcurves[nc].CriticalRadius(rmax,rmin,rave);
    std::cout << " critical curve radius: " << rmax/arcsecTOradians << " , " << rmin/arcsecTOradians
    << " , " << rave/arcsecTOradians << std::endl;
    
    //** make a Sersic source, there are a number of other ones that could be used
    SourceSersic source(22,0.1,0,1,0.5,z_source,zeropoint);
    source.setTheta(y[0].x);
    
    // adds the source image to the gridmap
    // The source redshift was fixed when the gridmap was made from the Lens.
    //  To change it use lens.ResetSourcePlane(zs) and remake the GridMap object.
    gridmap.AddSurfaceBrightnesses(&source);
    
    gridmap.writeFits(LensingVariable::SurfBrightness,"!surface_brightness.fits");
    
    // you can also make a Pixel map
    PixelMap pmap = gridmap.writePixelMap(LensingVariable::SurfBrightness);
    
    { // add another source
      SourceSersic source2(22,0.1,0,1,0.5,z_source,zeropoint);
      source2.setTheta(y[1].x);
      
      lens.ResetSourcePlane(5); // change redshift of source plane
      GridMap gridmap2(&lens,Ninit,center,range);
      
      gridmap2.AddSurfaceBrightnesses(&source2);
      pmap.AddGridMapBrightness(gridmap2);
    }

    pmap.printFITS("!surface_brightness2.fits");
   
    /*************************************************************
     Many other things are possible and easily done with GLAMER.
     Read the documentation for a more complete description of functionality.
     *************************************************************/
  }
 }
