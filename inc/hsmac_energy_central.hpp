
//  Df     df          df           df
// ---- = ---- + U * ------ + V * ------
//  Dt     dt         dksi         deta
//
// where:
//		U = a11 * u + a21 * v
//		V = a12 * u + a22 * v

         double FTX = (ap11(i,j)*ut + ap21(i,j)*vt) * tksi; 
         double FTY = (ap12(i,j)*ut + ap22(i,j)*vt) * teta; 

         // vt / J_t * (xksi * teta - xeta * tksi);
         // ur / J_t * (yeta * tksi - yksi * teta);
