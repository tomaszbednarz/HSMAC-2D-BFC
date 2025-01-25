
         double dif, o1;
 
         // higly optimalised UTOPIA scheme
         if ((i==1) || (i>=(NX-1))) {
            o1 = ( (uo(i+1,j)-uo(i-1,j)) - 1.0f*sgn(uo(i,j))*(uo(i+1,j)-2*uo(i,j)+uo(i-1,j)) ) / (2*dksi);
         } else {
            if (uo(i,j)>0.0f) dif=1.0f; else if (uo(i,j)<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * uo(i+2,j)
                  -(dif-2)/ 3 * uo(i+1,j)
                  +(dif  )/ 2 * uo(i  ,j)                     
                  -(dif+2)/ 3 * uo(i-1,j)
                  +(dif+1)/12 * uo(i-2,j) ) / dksi;
         }
         double FUX = uo(i,j) / J_u * (yeta * o1 - yksi * ueta);

         double vu = (vo(i,j)+vo(i+1,j)+vo(i,j-1)+vo(i+1,j-1)) / 4.0f;
         if ((j==1) || (j==NY)) {
            o1 = ( (uo(i,j+1)-uo(i,j-1)) - 1.0*sgn(vu)*(uo(i,j+1)-2*uo(i,j)+uo(i,j-1)) ) / (2*deta);
         } else {
            if (vu>0.0f) dif=1.0f; else if (vu<0.0f) dif=-1.0f; else dif=0.0f;
            o1= (  (dif-1)/12 * uo(i,j+2)
                  -(dif-2)/ 3 * uo(i,j+1)
                  +(dif  )/ 2 * uo(i,j  )                     
                  -(dif+2)/ 3 * uo(i,j-1)
                  +(dif+1)/12 * uo(i,j-2) ) / deta;
         }
         double FUY = vu / J_u * (xksi * o1 - xeta * uksi);
