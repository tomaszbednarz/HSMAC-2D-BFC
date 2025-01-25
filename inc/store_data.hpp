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


void read_data_array_double(ifstream *fp, Array<double,2> arr)
{
   Array<double,2>::const_iterator iter = arr.begin();
   Array<double,2>::const_iterator end  = arr.end();
   while (iter != end) {
      fp->read((char*)&(*iter),sizeof(double)); 
      ++iter;
   }
}

void store_data_array_double(ofstream *fp, Array<double,2> arr)
{
   Array<double,2>::const_iterator iter = arr.begin();
   Array<double,2>::const_iterator end  = arr.end();
   while (iter != end) {
      fp->write((char*)&(*iter),sizeof(double)); 
      ++iter;
   }
}

void load_stored_data(char *fn)
{
   ifstream fp (fn, ios::in | ios::binary);

   if ( !fp.is_open() )
   {
        cerr << "Unable to open to file: " << fn << endl;
        exit(1);
   }

   // skip header
   fp.seekg(25);
   
   char ver[50]; 
   fp.read((char*)ver,8); // = "2xxx0x0x";
   ver[8]=0; 
   cout << "*** reading flow data from file: " << fn << ", Version: " << ver << endl;

   fp.read((char*)&NX,sizeof(NX));
   fp.read((char*)&NY,sizeof(NY));

   cout << "\t- grid size: " << NX+1 << " x " << NY+1 << endl;

   // N-S equation parameters
   fp.read((char*)&AA,sizeof(AA));
   fp.read((char*)&BB,sizeof(BB));
   fp.read((char*)&CC,sizeof(CC));
   fp.read((char*)&DD,sizeof(DD));
   fp.read((char*)&T_0,sizeof(T_0));

   fp.read((char*)&iiter,sizeof(iiter));
   fp.read((char*)&niter,sizeof(niter));
   fp.read((char*)&iter,sizeof(iter));
   cout << "\t- starting iteration: " << iter << ", number of iterations computed: " << niter << endl;
   fp.read((char*)&dt,sizeof(dt));
   cout << "\t- delta time: " << dt << endl;

   fp.read((char*)&dksi,sizeof(dksi));
   fp.read((char*)&deta,sizeof(deta));
   cout << "\t- delta ksi: " << dksi << ", delta eta: " << deta << endl;

//   cout.precision(4);
//   cout.setf(ios::scientific,ios::floatfield);
   cout << "\t- AA: " << AA << ", BB: " << BB << ", CC: " << CC << ", DD: " << DD << ", T_0: " << T_0 << endl;

   init_tables(NX,NY);

   fp.read((char*)&pfactor,sizeof(pfactor));
   cout << "\t- pressure relaxation factor: " << pfactor << endl;

   cout << "\t- position x and position y" << endl;  
   read_data_array_double(&fp, x);
   read_data_array_double(&fp, y);
   cout << "\t- pressure" << endl;  
   read_data_array_double(&fp, p);
   cout << "\t- velocity -u and velocity -v components" << endl;  
   read_data_array_double(&fp, un);
   read_data_array_double(&fp, vn);
   cout << "\t- temperature" << endl;  
   read_data_array_double(&fp, tn);

/*
   read_data_array_double(&fp, uo);
   read_data_array_double(&fp, vo);
   read_data_array_double(&fp, to);

   cout << "\t- pressure-correction constants" << endl;  
   read_data_array_double(&fp, dppp);
   cout << "\t- Jacobians: center of the cells" << endl;  
   read_data_array_double(&fp, jacp);
   cout << "\t- Jacobians: -u velocity face" << endl;  
   read_data_array_double(&fp, jaceu);
   cout << "\t- Jacobians: -v velocity face" << endl;  
   read_data_array_double(&fp, jacnv);
*/
   init_pressure_velocity_correction();

   fp.close();
}

void store_data(char *fn)
{
   ofstream fp (fn, ios::out | ios::binary);

   if (fp.bad())
   {
        cerr << "Unable to write to file: " << fn << endl;
        exit(1);
   }

   cout << "*** writing flow data in file: " << fn << endl;

   // DO NOT CHANGE NEXT LINE!
   fp.write("HSMAC-BFC-2D-WARLOCK-VER-",25);
   // DO NOT CHANGE PREVIOUS LINE!

   // write version number (DATE)
   fp.write(VersionString,8);

   fp.write((char*)&NX,sizeof(NX));
   fp.write((char*)&NY,sizeof(NY));

   fp.write((char*)&AA,sizeof(AA));
   fp.write((char*)&BB,sizeof(BB));
   fp.write((char*)&CC,sizeof(CC));
   fp.write((char*)&DD,sizeof(DD));
   fp.write((char*)&T_0,sizeof(T_0));

   fp.write((char*)&iiter,sizeof(iiter));
   fp.write((char*)&niter,sizeof(niter));
   fp.write((char*)&iter,sizeof(iter));

   fp.write((char*)&dt,sizeof(dt));
   fp.write((char*)&dksi,sizeof(dksi));
   fp.write((char*)&deta,sizeof(deta));
   fp.write((char*)&pfactor,sizeof(pfactor));

   cout << "\t- position x and position y" << endl;  
   store_data_array_double(&fp, x);
   store_data_array_double(&fp, y);
   cout << "\t- pressure" << endl;  
   store_data_array_double(&fp, p);
   cout << "\t- velocity -u and velocity -v components" << endl;  
   store_data_array_double(&fp, un);
   store_data_array_double(&fp, vn);
   cout << "\t- temperature" << endl;  
   store_data_array_double(&fp, tn);
/*
   store_data_array_double(&fp, uo);
   store_data_array_double(&fp, vo);
   store_data_array_double(&fp, to);

   cout << "\t- pressure-correction constants" << endl;  
   store_data_array_double(&fp, dppp);
   cout << "\t- Jacobians: center of the cells" << endl;  
   store_data_array_double(&fp, jacp);
   cout << "\t- Jacobians: -u velocity face" << endl;  
   store_data_array_double(&fp, jaceu);
   cout << "\t- Jacobians: -v velocity face" << endl;  
   store_data_array_double(&fp, jacnv);
*/
   fp.close();
}


