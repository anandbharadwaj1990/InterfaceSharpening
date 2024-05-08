A simple FVM code to illustrate the algorithm of sharpening an interface in two-phase flows 
1. The fortran code needs to be compiled as shown: 
  gfortran -fdefault-real-8 2D_levelset_FVM_phi.f90

2. Run the code:
   ./a.out

Remarks 
1. LS.dat contains the interface marker points along with normals
2. 'Initial.csv' file is a smeared interface
3. 'Final.csv' file is the sharpened interface
4. '2d_phi.dat' file can be used to plot interface normals

  
