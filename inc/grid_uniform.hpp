//
// BFC + HSMAC method 
// 
// Programmed by Dr Tomasz BEDNARZ of JCU (School of Engineering)
// Email: tomasz.bednarz@jcu.edu.au
//
// November 2006, version 20061101.3
//
// All rights reserved (C)
//
// If you have this code, it means that it was given to you by
// the author.  If you use this code or any part of it for your
// own purposes, please mention the name of author and add his
// name to any papers or reports published with this ideas here.
// Also if you prefer F90, the code written in easy way to
// convert it to F90 - just copy and paste and change the 
// looping functions, etc...
// However, I recommend you to use C++ version due to many
// changes coming in the near future... 
//
// to compile the program you need to have:
// GNU C/C++
// Blitz++ library
//
// The code is solving NS equation on orthogonal meshes.
//
// REMEBER:
// The code is still under development process - please ask the
// me every 2nd month for the new and updated version.
//
// Things under developent:
// - corrections for non-orthogonality and highly skew grids
// - XML configuration files for universality
// - relative iteration progress notifications
// - particle transport 
// - postprocessing functions
//
//////////////////////////////////////////////////////////////////


void generate_uniform_grid(int sx, int sy)
{
   int i,j;
 
   NX = sx-1;
   NY = sy-1;

   init_tables(NX,NY);

double dup=0;
   for (j=0; j<=NY; j++) {
      for (i=0; i<=NX; i++) {
         x(i,j) = 0.03+LKSI*(double)(i)/NX+dup;
         y(i,j) = LETA*(double)(j)/NY;
//         cout << x << "  " << y << "\n";
      }
dup-=0.001;
   }

   // ----> compute points outside domain
   for (int i=0; i<=NX; i++) 
   {
      // bottom
      x(i,-1) = x(i,0) - ( x(i,1)-x(i,0) );
      y(i,-1) = y(i,0) - ( y(i,1)-y(i,0) );
      // top               
      x(i,NY+1) = x(i,NY) + ( x(i,NY)-x(i,NY-1) );
      y(i,NY+1) = y(i,NY) + ( y(i,NY)-y(i,NY-1) );
   }
   for (int j=-1; j<=(NY+1); j++) 
   {
      // left
      x(-1,j) = x(0,j) - ( x(1,j)-x(0,j) );
      y(-1,j) = y(0,j) - ( y(1,j)-y(0,j) );
      // right
      x(NX+1,j) = x(NX,j) + ( x(NX,j)-x(NX-1,j) );
      y(NX+1,j) = y(NX,j) + ( y(NX,j)-y(NX-1,j) );
   }
   // <----

   dksi = LKSI / NX;
//   dxi  = dksi;
   deta = LETA / NY;
//   det  = deta;

}
