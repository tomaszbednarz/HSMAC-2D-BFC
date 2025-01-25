
	double o1, dif;

	double v_at_u = (vo(i,j)+vo(i+1,j)+vo(i,j-1)+vo(i+1,j-1)) / 4.0f;
	double unew = ae11(i,j)*uo(i,j) + ae21(i,j)*v_at_u;
	double vnew = ae12(i,j)*uo(i,j) + ae22(i,j)*v_at_u;

	if ((i==1) || (i>=(NX-1))) {
	   o1 = ( (uo(i+1,j)-uo(i-1,j)) - 1.0f*sgn(unew)*(uo(i+1,j)-2*uo(i,j)+uo(i-1,j)) ) / (2*dksi);
        } 
        else 
        {
 	   if (unew>0.0f) dif=1.0f; else if (unew<0.0f) dif=-1.0f; else dif=0.0f;
           o1= (  (dif-1)/12 * uo(i+2,j)
                 -(dif-2)/ 3 * uo(i+1,j)
                 +(dif  )/ 2 * uo(i  ,j)                     
                 -(dif+2)/ 3 * uo(i-1,j)
                 +(dif+1)/12 * uo(i-2,j) ) / dksi;
	}
	double FUX = unew*o1;

	
	if ((j==1) || (j==NY)) {
	   o1 = ( (uo(i,j+1)-uo(i,j-1)) - 1.0*sgn(vnew)*(uo(i,j+1)-2*uo(i,j)+uo(i,j-1)) ) / (2*deta);
	} 
	else 
	{
	if (vnew>0.0f) dif=1.0f; else if (vnew<0.0f) dif=-1.0f; else dif=0.0f;
	   o1= (  (dif-1)/12 * uo(i,j+2)
                 -(dif-2)/ 3 * uo(i,j+1)
                 +(dif  )/ 2 * uo(i,j  )                     
                 -(dif+2)/ 3 * uo(i,j-1)
                 +(dif+1)/12 * uo(i,j-2) ) / deta;
	}
	double FUY = vnew*o1;
