

            double dif, o1;

            double uv = (uo(i,j)+uo(i,j+1)+uo(i-1,j)+uo(i-1,j+1)) / 4.0f;
            if ((i==1) || (i==NX)) {
                o1 = ( (vo(i+1,j)-vo(i-1,j)) - 1.0f*sgn(uv)*(vo(i+1,j)-2*vo(i,j)+vo(i-1,j)) ) / (2*dksi);
            } else {
                if (uv>0.0f) dif=1.0f; else if (uv<0.0f) dif=-1.0f; else dif=0.0f;
                o1= (  (dif-1)/12 * vo(i+2,j)
                      -(dif-2)/ 3 * vo(i+1,j)
                      +(dif  )/ 2 * vo(i  ,j)                     
                      -(dif+2)/ 3 * vo(i-1,j)
                      +(dif+1)/12 * vo(i-2,j) ) / dksi;
            }
            double FVX = uv / J_v * (yeta * o1 - yksi * veta);


            if ((j==1) || (j>=(NY-1))) {
                o1 = ( (vo(i,j+1)-vo(i,j-1)) - 1.0*sgn(vo(i,j))*(vo(i,j+1)-2*vo(i,j)+vo(i,j-1)) ) / (2*deta);
            } else {
                if (vo(i,j)>0.0f) dif=1.0f; else if (vo(i,j)<0.0f) dif=-1.0f; else dif=0.0f;
                o1= (  (dif-1)/12 * vo(i,j+2)
                      -(dif-2)/ 3 * vo(i,j+1)
                      +(dif  )/ 2 * vo(i,j  )                     
                      -(dif+2)/ 3 * vo(i,j-1)
                      +(dif+1)/12 * vo(i,j-2) ) / deta;
            }
            double FVY = vo(i,j) / J_v * (xksi * o1 - xeta * vksi);
