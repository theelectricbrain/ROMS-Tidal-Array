ROMS-Tidal-Array
================
###Warnings###
* This project is under development and subject to changes.
* It is a beta version and still suffers multiple bugs.
* It has been developed around an old version of ROMS (i.e. 3.6).


### Project description ###
#### Background:
* The tidal array representation method implemented in this modelling tool is based on a
  tidal turbine parameterisation method, accounting for the momentum capture as well as
  the sub-grid scale turbulence balance perturbations for each individual device.
  A detailed description and validation of this methodology can be found in Roc et al.
  2013[1] & 2014[2]. Each turbine behaviour is tuned via 4 parameters, namely Ct, Ctke
  and Cgls, respectively relating to the thrust coefficient, turbulent-cascade 
  short-circuiting and turbulence length-scale transfer[1].
* This tidal array representation method has been implemented within a large-scale ocean
  circulation modelling system that is the Regional Ocean Modeling System,
  [ROMS](www.myroms.org). ROMS is a member of a general class of three-dimensional,
  free-surface, terrain-following numerical models that solve the Reynolds-averaged
  Navier–Stokes equations using the hydrostatic and Boussinesq assumptions. In addition
  to the underlying hydro-dynamical engine – responsible for determining sea level height
  , and the three-dimensional circulation and transport of momentum, temperature and 
  salt – various turbulence closure sub-models are available for turbulent mixing. Note 
  that ROMS version 3.6 is used in the current distribution of ROMS-Tidal Array.

### Module Description ###
##### CPP Options
Three new CPP options have been added to the code and work as follows:
* TIDAL_TURBINE: switches on the tidal array computation. Note that the positions and
  features of the tidal turbines composing the array have to be specified in the
  “turbine input file”.
* I_ORIENTATION: makes all the tidal turbine rotors facing the I-direction of the grid
  (i.e. default orientation).
* J_ORIENTATION: makes all the tidal turbine rotors facing the J-direction of the grid
##### Input File
A new input file has been created in order to facilitate the use of this new
functionality. The input file (see “turbine_validation.in” in ROMS-Tidal-Array/cases/validation/) is composed of three main fields, “Lturbines”, “NTURBINES” and “POS”.
* “Lturbines” is a list of Boolean values (i.e. “T” or “F”) of “Ngrids” (i.e. number of nested and/or connected grid) element which switches the computation of turbine effects within nested and/or multiple connected grids. “T” would switch it on, “F” would switch it off.
* “NTURBINES” is a list of integer values of “Ngrids” (i.e. number of nested and/or connected grid) element which defines the number of turbines to account within nested and/or multiple connected grids.
* “POS” is composed of multiple lists, one of each turbine composing the array and one list per line. Each list/line is composed of 1 integer value (i.e. G) and 9 float values (i.e. xpos, ypos, zpos, Ct, Ctke, Cgls, Lc, Pa and Diam) defining the following turbines charcteristics:
  - “G”: Nested grid number
  - “xpos”: turbine x coordinate (in meters)
  - “ypos”: turbine y coordinate (in meters)
  - “zpos”: turbine z coordinate in depth ratio (i.e. hub height divided by water-column depth)
  - “Ct”: thrust coefficient
  - “Ctke”: Turbulent Kinetic Energy (TKE) correction parameter
  - “Cgls”: Generic Length Scale (GLS) correction parameter
  - “Lc”: blade chord length (m)
  - “Pa”: pitch angle (radians)
  - “Diam”: Turbine rotor diameter (m)

##### Output File
4 new fields have been added to the history output file which can be described as follows:
* “Dragfrc”: turbine induced momentum sink (m4.s-2)
* “Dragpwr”: turbine induced power sink (m5.s-3)
* “Dragtke”: turbine induced TKE correction (m2.s-2)
* “Draggls”: turbine induced GLS correction (m3.s-2)
Note that the tidal array module is compatible with the ROMS distributed-memory parallelism using OpenMPI libraries, besides the turbulence correction implementation is compatible with all the variations permitted by the GLS closure model (i.e. k-kl, k-ε, k-ω and k-generic).

### Version limitations ###
Version 0.9 of ROMS-Tidal Array has several known limitations:
* The tidal array functionality has been tested against all combination of options/couplings available in ROMS(please provide feedback, error log, ...)
* The current distribution is based on old version of ROMS and may therefore have possible compatibility problems with newer version of ROMS
* At this stage of development, modelled turbines work both ways with fixed directions and positions conditioned by the model grid.
* The formulations and coefficients of the turbine related terms are partially based on empirical tuning with publicly available database and therefore may require custom modifications as well as further improvements as validation data becomes more and more available.

### Best practises ###
The following recommendations will help enjoying the ROMS-Tidal Array experience to the fullest. The learning curve is steep but worthwhile.
* First of all, join the [ROMS/TOMS group](https://www.myroms.org/index.php?page=login). Their community forum will quickly become your best ally when inevitably battling with unknown error messages. Do not forget to acknowledge their amazing work when [disseminating your work](https://www.myroms.org/index.php?page=License_ROMS) and try as much as possible to actively participate to the community effort.
* Successively set-up, run and validate the “upwelling” test case provide by ROMS web site before to try out the “TIDAL_TURBINE” functionality.
* Successively run test case provided in “~/ROMS-Tidal Array/cases/validation/” and compare the so-generated history output file with “ocean_validation_his.nc”
* Once you have reached this point, you should be able to develop your own scenarios yet remember to make (x, y) turbine hub coordinates with rho-point coordinates and note that experience has shown that the horizontal spatial resolution should not exceed a 1/3 of the diameter along the turbine rotor and 1 diameter across the turbine rotor1.

### References ###
1. Roc et al. 2013. "Methodology for tidal turbine representation in ocean circulation". Renewable Energy, Volume 51, pg. 448-464
2. Roc et al. 2014: "Tidal turbine representation in an ocean circulation model: Towards realistc applications". Ocean Engineering, Volume 78, pg. 95-111

### Acknowledgement ###
* ROMS/TOMS group  for their amazing work on the ROMS model (see www.myroms.org)
* IT Power Consulting Ltd for supporting and financing this project (see www.itpower.co.uk)
* Electric Brain for administrating and mainting this repository/group ( see electricbrain.fr)
* Special acknowledgement for the contributions of [Alice Jane Goward Brown](https://www.bangor.ac.uk/oceansciences/staff/phd-students/alice-goward-brown) and [Elizabeth Brasseale](http://www.ocean.washington.edu/home/Elizabeth+Brasseale)

### Legal Information ###
* ROMS core code: Copyright (c) 2002-2016 The ROMS/TOMS Group, MIT/X License (See http://www.opensource.org/licenses/mit-license.php)
* ROMS tidal array module: Copyright (c) 2016 Thomas Roc and ITPower, Licensed under an Affero GPL style license v3.0 (see License_ROMS-Tidal Array.txt)
