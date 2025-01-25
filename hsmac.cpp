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
// the author and you are not allowed to give it to others.  
// If you use this code or any part of it for your own purposes, 
// please mention the name of author and add his
// name to any papers or reports published with this ideas here.
//
// Also if you prefer F90, the code written in easy way to
// convert it to F90 - just copy and paste and change the 
// looping functions, etc...
//
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

#include <blitz/array.h>   

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace blitz;
using namespace std;

#include "./inc/global.hpp"
#include "./inc/hsmac_energy.hpp"
#include "./inc/hsmac_momentum.hpp"
#include "./inc/hsmac_pressure_second.hpp"
#include "./inc/tecplot.hpp"
#include "./inc/store_data.hpp"
#include "./inc/grid_uniform.hpp"
#include "./inc/bc_temperature.hpp"

void init_tables(int nx, int ny)
{
   x.resize(Range(-1,nx+1),Range(-1,ny+1));
   y.resize(Range(-1,nx+1),Range(-1,ny+1));

   uo.resize(Range(0,nx  ),Range(0,ny+1));
   un.resize(Range(0,nx  ),Range(0,ny+1));

   vo.resize(Range(0,nx+1),Range(0,ny  ));
   vn.resize(Range(0,nx+1),Range(0,ny  ));

   to.resize (Range(0,nx+1),Range(0,ny+1));
   tn.resize (Range(0,nx+1),Range(0,ny+1));

   p.resize (Range(0,nx+1),Range(0,ny+1));

   psi.resize (Range(0,nx+1),Range(0,ny+1));

   dppp.resize (Range(0,nx+1),Range(0,ny+1));

   jacp.resize  (Range(0,nx+1),Range(0,ny+1));
   jaceu.resize (Range(0,nx+1),Range(0,ny+1));
   jacnv.resize (Range(0,nx+1),Range(0,ny+1));

   ap11.resize (Range(0,nx+1),Range(0,ny+1));
   ap12.resize (Range(0,nx+1),Range(0,ny+1));
   ap21.resize (Range(0,nx+1),Range(0,ny+1));
   ap22.resize (Range(0,nx+1),Range(0,ny+1));

   ae11.resize (Range(0,nx+1),Range(0,ny+1));
   ae12.resize (Range(0,nx+1),Range(0,ny+1));
   ae21.resize (Range(0,nx+1),Range(0,ny+1));
   ae22.resize (Range(0,nx+1),Range(0,ny+1));

   an11.resize (Range(0,nx+1),Range(0,ny+1));
   an12.resize (Range(0,nx+1),Range(0,ny+1));
   an21.resize (Range(0,nx+1),Range(0,ny+1));
   an22.resize (Range(0,nx+1),Range(0,ny+1));

   // for 2nd derivatives
   ap2d1.resize (Range(1,nx),Range(1,ny));
   ap2d2.resize (Range(1,nx),Range(1,ny));
   ae2d1.resize (Range(1,nx),Range(1,ny));
   ae2d2.resize (Range(1,nx),Range(1,ny));
   an2d1.resize (Range(1,nx),Range(1,ny));
   an2d2.resize (Range(1,nx),Range(1,ny));

   ap11 = 0.0f;
   ap12 = 0.0f;
   ap21 = 0.0f;
   ap22 = 0.0f;

   ae11 = 0.0f;
   ae12 = 0.0f;
   ae21 = 0.0f;
   ae22 = 0.0f;

   an11 = 0.0f;
   an12 = 0.0f;
   an21 = 0.0f;
   an22 = 0.0f;


   jacp = 0.0f;
   jaceu = 0.0f;
   jacnv = 0.0f;
   
   dppp = 0.0f;
   p    = 0.0f;
   uo   = 0.0f;
   vo   = 0.0f;
   un   = 0.0f;
   vn   = 0.0f;
   tn   = 290.0f;
   to   = 290.0f;
   psi  = 0.0f;
}



void read_grid_xy(char *f1, char *f2)
{
  
   ifstream fx (f1);
   ifstream fy (f2);

   if ( fx.is_open() ) {
      fx >> NX;
      NX--;
   }
   if ( fy.is_open() ) {
      fy >> NY;
      NY--;
   }

   init_tables(NX,NY);

   for (int j=0; j<=NY; j++) 
   {
      double yyy;
      fy >> yyy;
      // cout << yyy << "\n";
      for (int i=0; i<=NX; i++) 
      {
        y(i,j) = yyy;
      }
   }

   for (int i=0; i<=NX; i++) 
   {
      double xxx;
      fx >> xxx;
      // cout << yyy << "\n";
      for (int j=0; j<=NY; j++) 
      {
        x(i,j) = xxx;
      }
   }

   fx.close();
   fy.close();

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


void read_grid(char *filename)
{
   std::istringstream sx,sy,stm;
   string line;
   string s1;
   double d1, d2;

   ifstream fp (filename);

   if ( fp.is_open() )
   {
      getline (fp,line);
      getline (fp,line);

      // parse X size
      s1 = line.erase(0,line.find("=",0)+1);
      sx.str( line.substr(0,line.find(",",0)) );
      sx >> NX;
      // parse Y size
      s1 = line.erase(0,line.find("=",0)+1);
      sy.str( s1 );
      sy >> NY;

      NX--; // denotes number of cells as well
      NY--;
   
      cout << "read grid (no of cells) -> NX = " << NX << ", NY = " << NY << "\n";

      init_tables(NX,NY);
   
      for (int j=0; j<=NY; j++) 
      {
         for (int i=0; i<=NX; i++) 
         {
            getline(fp,line);
            stm.str(line);
            stm >> x(i,j) >> y(i,j);
            x(i,j) = x(i,j)*0.03; //66.66666666666666;
            y(i,j) = y(i,j)*0.03; //66.66666666666666;
         }
      }

      fp.close();
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



void boundary_velocity()
{
   int i,j;

   // ---------------------------> u
   for (j=1; j<=NY; j++) 
   {
      un(NX,j) = 0.0f; // left
      un(0, j) = 0.0f; // right
   }
   for (i=0; i<=NX; i++)
   {
//      un(i,NY+1) = -1.0f; ///1000.0f; // top

///////      un(i,NY+1) = 2* ( 1.0f ) - un(i,NY); ///1000.0f; // top
/*
if (iter>=1000)
{
      un(i,NY+1) = 1.0f; ///1000.0f; // top
}
if (iter<1000)
{
      un(i,NY+1) = 1.0f/1000.0f*(double)iter; ///1000.0f; // top
}
*/
      un(i,0   ) = -un(i,1);      // bottom

//un(i,0)=0.0f;
un(i,NY+1)=-un(i,NY);
   }
   // ---------------------------> v
   for (j=1; j<NY; j++)
   {
      vn(NX+1,j) = -vn(NX,j); // right
      vn(0   ,j) = -vn(1,j ); // left
//vn(NX+1,j)=0.0f;
//vn(0,j)=0.0f;
   }
   for (i=0; i<=(NX+1); i++)
   {
      vn(i,NY) = 0.0f; // top
      vn(i,0 ) = 0.0f; // bottom
//vn(i,NY)=0.0f;
//vn(i,0)=0.0f;
   }

/*
   // U left and right
   for (j=1; j<jmax; j++) {
      //un[im1][j]=un[im2][j];               //0.0f;
      //un[0][j]=1.0f;                       //0.0f;
      un[im1][j]=0.0f;
      un[0][j]=0.0f;
   }
   // U top and bottom
   for (i=0; i<imax; i++) {
      un[i][jmax]=-un[i][jm1];
      un[i][0]=-un[i][1];
   }
   // V left and right
   for (j=0; j<jm1; j++) {
      vn[imax][j]=-vn[im1][j];
      vn[0][j]=-vn[1][j];
   }
   // V top and bottom
   for (i=0; i<=imax; i++) {
      vn[i][jm1]=0.0f;
      vn[i][0]=0.0f;
   }

*/

}

void boundary_pressure()
{
   int i,j;

   for (j=1; j<=NY; j++)
   {
      p(0   ,j) = p(1 ,j); // left
      p(NX+1,j) = p(NX,j); // right
   }
   for (i=0; i<=(NX+1); i++)
   {
      p(i,NY+1) = p(i,NY); // top
      p(i,0   ) = p(i,1 ); // bottom
   }
}



void calculate_psi()
{
   int i,j;

   for (i=1; i<=NX; i++) psi(i,0) = psi(i-1,0) - ( (vn(i,0)+vn(i+1,0))/2.0f ) * dksi; 

   for (i=0; i<=NX; i++)
   {
      for (j=1; j<=NY; j++)
      {
         psi(i,j) = psi(i,j-1) + ( (un(i,j)+un(i,j+1))/2.0f ) * deta;
      }
   }
}


double calculate_eps(Array<double,2> qn, Array<double,2> qo)
{
   double eps;

   eps = 0.0f;
 
   for (int j=1; j<NY; j++) {
      for (int i=1; i<NX; i++) {
         double e = fabs( qn(i,j)-qo(i,j) );
         if (e > eps) {
            eps = e;
         }
      }
   }
   return eps;
}

void read_config()
{
   string line;   
   string f1,f2;
   
   int is_new_file;
   int tx,ty;

   ifstream fp ("input\\hsmac.in", ios::in | ios::binary);

   if ( fp.is_open() )
   {
      // new case? 0-new, 1-old
      getline (fp,line);
      fp >> is_new_file;

      if (is_new_file) {

         getline (fp,line);
         getline (fp,line);
         getline (fp,line);
         getline (fp,line);
         getline (fp,line);

         int grid;
         fp >> grid;
        
         // new case: grid 0-one file, 1-two files (x,y), 2-uniform grid 
         switch (grid) {
            case 0:
               fp >> f1;
               read_grid((char*)f1.c_str());
               break;
            case 1:
               fp >> f1 >> f2;
               read_grid_xy((char*)f1.c_str(),(char*)f2.c_str());
               break;
            case 2:
               fp >> tx >> ty;
               generate_uniform_grid(tx,ty);
               break;
         }

         getline (fp,line);
         getline (fp,line);
         fp >> niter >> iiter;
         cout << "Number of iterations: " << niter << endl;
         cout << "Number of inner iterations (pressure): " << iiter << endl;

         getline (fp,line);
         getline (fp,line);
         fp >> pfactor;
         cout << "Relaxation factor for pressure: " << pfactor << endl;

         getline (fp,line);
         getline (fp,line);
         getline (fp,line);
         getline (fp,line);
         getline (fp,line);
         getline (fp,line);

         // dimensional: AA=1/rho, BB=mu/rho, CC=g*beta, DD=alpha=k/(rho*cp), T_0=temp_ref
         fp >> AA >> BB >> CC >> DD >> T_0;
         cout << "AA: " << AA << ", BB: " << BB << ", CC: " << CC << ", DD: " << DD << ", T_0: " << T_0 << endl;

         getline (fp,line);
         getline (fp,line);
         fp >> dt;
         cout << "dt: " << dt << endl;
        
         iter = 0;

         boundary_velocity();
         boundary_temperature();

         init_pressure_velocity_correction();

//         init_pressure_velocity_correction();
//         init_pressure_velocity_correction_dwa();

      } else {
         // read existing file (grid, previous time steps, etc.)
         getline (fp,line);
         getline (fp,line);
         fp >> line;
         load_stored_data((char*)line.c_str());
      }

      fp.close(); 
   } else {
      cerr << "Unable to open config file...\n" << endl;
      exit(1);
   }
}

int main()
{
   read_config();

   int l=0,k=0;

   do {
      // advance
      cycleArrays(uo, un);
      cycleArrays(vo, vn);
      cycleArrays(to, tn);

      calculate_velocity();

      pressure_velocity_correction();

      calculate_energy();

      iter++;

      if ( (iter % 100) == 0 )
      {
         cout.precision(4);
         cout.setf(ios::scientific,ios::floatfield);

         if (l==0)
         cout << setw(8) << "iter" << "  continuity" << "  x-velocity" << "  y-velocity" << "  temperatur"<< "  innerp" << "  time" << endl; // << "  "<< abs(dmax) << "\t" << inner << "\t" << iter*dt << "\n";
         cout << setw(8) << iter << "  " << abs(dmax) << "  " << calculate_eps(un,uo) << "  " << calculate_eps(vn,vo) << "  " << calculate_eps(tn,to) << "  " << setw(6) << inner << "  " << fixed << iter*dt << endl;
   
         l++;
         if (l>10) l=0;
      }

      // tecplot
      if ( ((iter % 250) == 0) && (iter !=0) ) //20000
      {
         char buffer [50];
         char buffer2[50];
         sprintf (buffer , "results\\iter%09d.plt",iter);
         sprintf (buffer2, "time=%.5f",iter*dt);
         //cout << buffer << "\n";
         calculate_psi();
         write_tecplot(buffer, "ufff", buffer2);
      }


   } while (++k < niter);


   calculate_psi();
   write_tecplot("plik.plt", "tomcio", "ZONE 001");

   store_data("results\\_20000.wlk");

/*
   ofstream fp ("gridaaaaa.dat");
   if ( fp.is_open() )
   {
      fp << "VARIABLES = x, y, u, v, t, psi\n";
      fp << "ZONE I=" << NX+1 << ", J=" << NY+1 << " F=POINT\n";
      fp.setf(ios::showpoint | ios::right);
      fp.precision(15);
      for (int j=0; j<=(NY); j++) {
         for (int i=0; i<=(NX); i++) {
            fp.setf( ios::scientific);
            fp.precision(8);
            fp << "  " << x(i,j) << "  " << y(i,j) << "  " << (un(i,j)+un(i,j+1))/2.0f << "  " << (vn(i,j)+vn(i+1,j))/2.0f << "  " << ((tn(i,j)+tn(i+1,j)+tn(i,j+1)+tn(i+1,j+1))/4.0f) << "  " << psi(i,j) << "\n";          
//            fp << "  " << x(i,j) << "  " << y(i,j) << "  " << un(i,j) << "  " << vn(i,j) << "  " << tn(i,j) << "  " << psi(i,j) << "\n";          
//            //fprintf (fp,"%-3.16f  %-3.16f\n",x(i,j),y(i,j));
//            fp << "  " << x(i,j) << "  " << y(i,j) << endl;          
         }
      }
      fp.close();
   }
*/

   return 0;
}

