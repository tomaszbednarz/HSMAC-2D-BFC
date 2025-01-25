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


void write_tecplot_name(ofstream *fp, char *name)
{
   int i, tmp;

   i=0;
   do {
      tmp = name[i];
      fp->write((char*)&tmp,sizeof(tmp));
      i++;
   } while(name[i]!=0);
   tmp = 0;
   fp->write((char*)&tmp,sizeof(tmp));
}

void write_tecplot(char* plik, char* title, char* zone_name)
{
   char *nag = "#!TDV102";
   int tmp = 1;
   char *var_x = "x";
   char *var_y = "y";
   char *var_u = "u";
   char *var_v = "v";
   char *var_t = "t";
   char *var_p = "p";
   char *var_psi = "psi";
   float zone_marker = 299.0f;
   float ftmp;

   ofstream fp (plik, ios::out | ios::binary);

   if ( fp.is_open() )
   {
      // i. magic number, version number
      fp.write (nag,8);

      // ii. integer value of 1 (byte order)
      fp.write((char*)&tmp,sizeof(tmp));
      
      // iii. the TITLE and variable names
      write_tecplot_name(&fp, title);

      // number of variables
      tmp = 7;
      fp.write((char*)&tmp,sizeof(tmp));
      // variable names
      write_tecplot_name(&fp, var_x);
      write_tecplot_name(&fp, var_y);
      write_tecplot_name(&fp, var_u);
      write_tecplot_name(&fp, var_v);
      write_tecplot_name(&fp, var_t);
      write_tecplot_name(&fp, var_p);
      write_tecplot_name(&fp, var_psi);

      // iv. ZONES 
      // zone marker
      fp.write((char*)&zone_marker,sizeof(zone_marker));
      // zone name
      write_tecplot_name(&fp, zone_name);
      // zone color
      tmp = -1;
      fp.write((char*)&tmp,sizeof(tmp));
      // ZoneType 0=ORDERED, 1=FELINESEG, ...
      tmp = 0;
      fp.write((char*)&tmp,sizeof(tmp));
      // DataPacking 0=Block, 1=Point
      tmp = 1;
      fp.write((char*)&tmp,sizeof(tmp));
      // Specify Var Location. 0 = Don't specify, all date is located at the nodes, 1 = ...
      tmp = 0;
      fp.write((char*)&tmp,sizeof(tmp));
      // Number of user defined face neighbor connections (vale >= 0)
      tmp = 0;
      fp.write((char*)&tmp,sizeof(tmp));
      // IMax
      tmp = NX+1;
      fp.write((char*)&tmp,sizeof(tmp));
      // JMax
      tmp = NY+1;
      fp.write((char*)&tmp,sizeof(tmp));
      // Auxiliary data name/value pair
      tmp = 1;
      fp.write((char*)&tmp,sizeof(tmp));
      tmp = 0;
      fp.write((char*)&tmp,sizeof(tmp));

      // v. Geometries
      // skip geometry defs
      ftmp = 357.0f;
      fp.write((char*)&ftmp,sizeof(ftmp));

      // data section
      ftmp = 299.0f;
      fp.write((char*)&ftmp,sizeof(ftmp));
      // var data format 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
      tmp = 1;
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      fp.write((char*)&tmp,sizeof(tmp));
      // has variable sharing 0=no, 1=yes
      tmp = 0;      
      fp.write((char*)&tmp,sizeof(tmp));
     
      // Zone number to share connectivity list with (-1 = no sharing)
      tmp = -1;
      fp.write((char*)&tmp,sizeof(tmp));
      
      for (int j=0; j<=(NY); j++) {
         for (int i=0; i<=(NX); i++) {
            ftmp = x(i,j);
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = y(i,j);
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = (un(i,j)+un(i,j+1))/2.0f;
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = (vn(i,j)+vn(i+1,j))/2.0f;
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = ((tn(i,j)+tn(i+1,j)+tn(i,j+1)+tn(i+1,j+1))/4.0f);
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = ((p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))/4.0f);
            fp.write((char*)&ftmp,sizeof(ftmp));
            ftmp = psi(i,j);
            fp.write((char*)&ftmp,sizeof(ftmp));
         }
      }

      fp.close();
   }
}
