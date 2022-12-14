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

  /********************
   set parameter file name:
   A default name is used or a name is taken as a command line argument
   *********************/
/*
  std::string paramfile;
  std::cout << "initializing model" << std::endl;
  //string paramfile;
  if(arg > 1) paramfile.assign(argv[1],strlen(argv[1]));
  else paramfile = "sample_paramfile";
  std::cout << "using parameter file: " << paramfile << std::endl;
  */
  // read parameter file
  //InputParams params(paramfile);

  
  //****************************************
  // **** random number seed, not actually used in this example, but
  // ***** needed for lens constructor
  //****************************************
  long seed = -1827675;

  //**** set the cosmology  ******
  COSMOLOGY cosmo(CosmoParamSet::Planck18);

  
  // construct some massive objects call "halos"
  
  
  // an NFW halo - mass,Rmax,redshift,concentration,axis ratio
  //               ,position angle,number of stars
  
  double z_lens = 0.3; // redshift of lens
  //LensHaloNFW halo_nfw(1.0e13,1.0,z_lens,4,1,0,0);
  LensHaloNFW halo_nfw(1.0e13,1.0,z_lens,4,1,0,cosmo);
  halo_nfw.setTheta(0,0); // set angular position
  // a NonSingular Isothermal Ellipsoid - mass,redshift,sigma
  // , core radius,axis ratio,position angle,number of stars
  LensHaloRealNSIE halo_nsie(8.0e11,z_lens, 300, 0, 0.5, 0, cosmo);
  halo_nsie.setTheta(0,0);

  // construct lens
  Lens lens(&seed,2.0,cosmo);
  
  // insert object into the lens
  lens.insertMainHalo(halo_nfw,true);
  lens.insertMainHalo(halo_nsie,true);
  
  std::cout << "Lens constructed" << std::endl;
  
  /*******************************
  
   Now that we
   
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
  {
    std::cout << "Constructing initial GridMap ..." << std::endl;
    GridMap gridmap(&lens,Ninit,center,range);
    std::cout << "constructed" << std::endl;

    /***
     make an images of kappa
     The center,range and pixel numbers do not have to match the grid, but 
     in this case they are set to match.
     ****/
    gridmap.writeFitsUniform(center,gridmap.getInitNgrid()
                           ,gridmap.getInitNgrid(),LensingVariable::KAPPA,"!gridmap");

  }
  /****************************
   Here a Grid is constructed which, unlike GridMap can be refined for finding images, etc.
   ****************************/
  
  std::cout << "Constructing initial Grid ..." << std::endl;
  Grid grid(&lens,Ninit,center,range/2);
  std::cout << "constructed" << std::endl;
  
  /*
   There we make same maps of the lensing quantities.
   This could be done with GridMap since no refinement has been done yet.
   The "!" in front of the file name causes it to overwrite a file with that name.  Suffixes are added (e.g. .kappa.fits).
   */
  
  grid.writeFits(2,LensingVariable::KAPPA,"!initgrid");
  grid.writeFits(2,LensingVariable::ALPHA,"!initgrid");
  grid.writeFits(2,LensingVariable::GAMMA,"!initgrid");
  grid.writeFits(2,LensingVariable::INVMAG,"!initgrid");
  
  PixelMap map = grid.writePixelMap(center,grid.getInitNgrid(), grid.getInitRange()/grid.getInitNgrid() ,LensingVariable::ALPHA);

  /******************************************
   Now we are going to look for same caustics
   ******************************************/

  // The CriticalCurve class contains information about a caustic
  std::vector<ImageFinding::CriticalCurve> critcurves(100);
  int Ncrit;  // number of caustics
  // resolution to which the critical curves should be refined (radians)
  PosType resolution = 0.1*arcsecTOradians;
  
  std::cout << "Looking for critical curves ..." << std::endl;
  // Find all the critical curves
  ImageFinding::find_crit(&lens,&grid,critcurves, &Ncrit,resolution,0.0,true);
  std::cout << Ncrit << " critical curves found." << std::endl;

  PosType Xrange[2]={0,0},Yrange[2]={0,0};
  PosType rmax,rmin,rave;
  for(auto it : critcurves){
    it.CriticalRadius(rmax,rmin,rave);
    std::cout << (int)(it.type) << " " << it.critical_area << " " << it.caustic_area << std::endl;;
  }
  if(Ncrit > 0){
    Point_2d p1,p2;
    
    // find a box on the image plane that contains all of the critical curves
    for(int ii = 0;ii< Ncrit;++ii){
      critcurves[ii].CritRange(p1,p2);
 
      Xrange[0] = MIN(Xrange[0],p1[0]);
      Xrange[1] = MAX(Xrange[1],p2[0]);
      
      Yrange[0] = MIN(Yrange[0],p1[1]);
      Yrange[1] = MAX(Yrange[1],p2[1]);
    }
    
    // make a PixelMap which is used for IO of images
    PixelMap map(center,2*grid.getInitNgrid(),grid.getInitRange()/grid.getInitNgrid()/2);
    
    // this draws the critical curves on the image
    for(int ii = 0;ii< Ncrit;++ii){
      map.AddCurve(critcurves[ii].critcurve,ii+1);
      map.AddCurve(critcurves[ii].caustic_curve_intersecting,ii+1);
    }

    // print an image of the caustic
    map.printFITS("!caustics_crit.fits");
    
    /**************************
     Now we are going to pick a caustic, put a source 
     in it and find its images
     **************************/
  
    std::vector<Point_2d> ys;  // source positions
    // this is a random number generator
    Utilities::RandomNumbers_NR rng(seed);
    
    // This finds random points within this caustic.
    // In this case it is just one point.
    int nc=0;  // pick a caustic
    for (int i=1 ; i < critcurves.size() ; ++i) {
      if(critcurves[nc].critical_area > critcurves[i].critical_area) nc = i;
    }
    critcurves[nc].RandomSourceWithinCaustic(1,ys,rng);

    // ImageInfo is a class that contains information about images
    std::vector<ImageInfo> imageinfo;
    int Nimages;  // number of images that will be found
    size_t Nimagepoints;  // total number of points within the images
    
    // This finds the images made by a source at ys[0] and stores
    // information about them in imageinfo.
    // These are just circular sources with constant surface brightness
    ImageFinding::find_images_kist(&lens,ys[0].x,0.05*arcsecTOradians
                    ,&grid,&Nimages,imageinfo,&Nimagepoints,
                            0, true, 2);
    
    std::cout << "Number of images: " << Nimages << std::endl;
    
    // print an image of the caustic
    map.printFITS("!test.caustics.fits");
    
    map.Clean();
    // add images to the PixelMap
    map.AddImages(imageinfo.data(),Nimages,0);
    
    // output the PixelMap as a fits file
    map.printFITS("!test.image.fits");
  }
  
  // If the grid is now output at twice the original resolution
  // some additional structure might be see because of the refinement
  grid.writeFits(center,2*grid.getInitNgrid(),grid.getInitRange()/grid.getInitNgrid()/2, LensingVariable::INVMAG,"!initgrid_refined");
  grid.writeFits(center,2*grid.getInitNgrid(),grid.getInitRange()/grid.getInitNgrid()/2, LensingVariable::KAPPA,"!initgrid_refined");
  
  
  //*************************************************************
  //**** Let put a a more complicated source in the image
  //****************************************
  
  // select the largest area tangential critical curve
  int i=0;
  for (int ii=1 ; ii < critcurves.size() ; ++ii) {
    if(critcurves[ii].critical_area > critcurves[i].critical_area) i = ii;
  }

  //while(critcurves[i].type != tangential) ++i;

  //*** find a source position within the tangential caustic
  Utilities::RandomNumbers_NR random(seed);   //*** random number generator
  std::vector<Point_2d> y;                    //*** vector for source positions
  critcurves[i].RandomSourceWithinCaustic(1,y,random); //*** get random points within first caustic
  
  critcurves[i].CriticalRadius(rmax,rmin,rave);
  std::cout << " critical curve radius: " << rmax/arcsecTOradians << " , " << rmin/arcsecTOradians
  << " , " << rave/arcsecTOradians << std::endl;
  
  PosType zs = 2; //** redshift of source
  //** make a Sersic source, there are a number of other ones that could be used
  SourceSersic source(22,0.1,0,1,0.5,zs);
  source.setTheta(y[0].x);
  
  /** reset the source plane in the lens from the one given in the
   parameter file to this source's redshift
   */
  lens.ResetSourcePlane(zs,false);
  std::vector<ImageInfo> imageinfo;
  int Nimages;
  Point_2d p1,p2;
  critcurves[i].CritRange(p1,p2);  // this is to make the image fit the critical curve


  // get information on observational / telescope parameters

  Observation ob(Telescope::KiDS_i,100, 100);     // KiDS-like observation
 
  // high res image
  PixelMap mapI(critcurves[i].critical_center.x,200,ob.getPixelSize()/3);
  std::cout << "Mapping source ..." << std::endl;
  
  /*** there are two different ways to map the images
   ImageFinding::map_images() will adapt the grid to make the image smooth
   ImageFinding::map_images_fixedgrid() will use the grid as is without further ray shooting
   */
  ImageFinding::map_images(&lens,&source,&grid,&Nimages,imageinfo,source.getRadius()
                           ,source.getRadius()/100,0,ExitCriterion::EachImage,false,true);
  
  //ImageFinding::map_images_fixedgrid(&source,&grid,&Nimages,imageinfo
  //,source.getRadius(),true,true);
  
  //*** add sources images to the plot.
  mapI.AddImages(imageinfo,Nimages);
  mapI.printFITS("!image.fits");  // output image without psf and noise

  
  // new map with resolution of observation
  PixelMap mapII(critcurves[i].critical_center.x,100,ob.getPixelSize());
  
   // copy other image in taking into account the different sized pixels
  mapII.copy_in(mapI);

  // now add a unlensed source as the lens galaxy
  SourceSersic lens_galaxy(22,0.2,45*degreesTOradians,4,0.5,z_lens);
  
  mapII.AddSource(lens_galaxy);  // add source to image
  
  Utilities::RandomNumbers_NR ran(seed);
  
  ob.Convert(mapII,true, true, ran);  // apply psf and noise

  mapII.printFITS("!image_withnoise.fits");

 
  
  /*************************************************************
   Many other things are possible and easily done with GLAMER.
   Read the documentation for a more complete description of functionality.
   *************************************************************/

 }
