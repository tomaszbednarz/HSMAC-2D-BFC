

	double dif, o1;

	double u_at_v = (uo(i,j)+uo(i,j+1)+uo(i-1,j)+uo(i-1,j+1)) / 4.0f;
	
	double unew = an11(i,j)*u_at_v + an21(i,j)*vo(i,j);
	double vnew = an12(i,j)*u_at_v + an22(i,j)*vo(i,j);

	if ((i==1) || (i==NX)) {
           o1 = ( (vo(i+1,j)-vo(i-1,j)) - 1.0f*sgn(unew)*(vo(i+1,j)-2*vo(i,j)+vo(i-1,j)) ) / (2*dksi);
        } else {
           if (unew>0.0f) dif=1.0f; else if (unew<0.0f) dif=-1.0f; else dif=0.0f;
              o1= (  (dif-1)/12 * vo(i+2,j)
                    -(dif-2)/ 3 * vo(i+1,j)
                    +(dif  )/ 2 * vo(i  ,j)                     
                    -(dif+2)/ 3 * vo(i-1,j)
                    +(dif+1)/12 * vo(i-2,j) ) / dksi;
        }

	double FVX = unew*o1;

	
	if ((j==1) || (j>=(NY-1))) {
	   o1 = ( (vo(i,j+1)-vo(i,j-1)) - 1.0*sgn(vnew)*(vo(i,j+1)-2*vo(i,j)+vo(i,j-1)) ) / (2*deta);
	} else {
	   if (vnew>0.0f) dif=1.0f; else if (vnew<0.0f) dif=-1.0f; else dif=0.0f;
	      o1= (  (dif-1)/12 * vo(i,j+2)
                    -(dif-2)/ 3 * vo(i,j+1)
                    +(dif  )/ 2 * vo(i,j  )                     
                    -(dif+2)/ 3 * vo(i,j-1)
                    +(dif+1)/12 * vo(i,j-2) ) / deta;
        }
	double FVY = vnew*o1;

