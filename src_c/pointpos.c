/****************************************************************************

Routine:	       	pointpos.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	July 1988, June 1989
	
Description:	       	point positioning with least squares adjustment
		       	(multimedia, 2 or 3 channels)
	
Routines contained:		

****************************************************************************/
#include "ptv.h"
#include "lsqadj.h"

void dist_to_ray(x, y, Ex, I, G, ap, mm, Xp,Yp,Zp, dist)

Exterior	Ex;
Interior	I;
Glass   	G;
ap_52		ap;
mm_np		mm;
double		x, y,Xp,Yp,Zp, *dist;
{
  double    gX[4],gY[4],gZ[4],a[4],b[4],c[4];
  double    x01,x02,x03,x12,x13,x23;
  double    y01,y02,y03,y12,y13,y23;
  double    z01,z02,z03,z12,z13,z23;
  
    
  	
  *dist=1;

}
 
void pos_from_ray(Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp, dist)

Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp, *dist;
{
  double    gX[4],gY[4],gZ[4],a[4],b[4],c[4];
  double    x01,x02,x03,x12,x13,x23;
  double    y01,y02,y03,y12,y13,y23;
  double    z01,z02,z03,z12,z13,z23;
  
    
  ray_tracing_v2 (x1, y1, Ex[0], I[0], G[0], mmp, &gX[0], &gY[0], &gZ[0], &a[0], &b[0], &c[0]);
  ray_tracing_v2 (x2, y2, Ex[1], I[1], G[1], mmp, &gX[1], &gY[1], &gZ[1], &a[1], &b[1], &c[1]);
  ray_tracing_v2 (x3, y3, Ex[2], I[2], G[2], mmp, &gX[2], &gY[2], &gZ[2], &a[2], &b[2], &c[2]);
  ray_tracing_v2 (x4, y4, Ex[3], I[3], G[3], mmp, &gX[3], &gY[3], &gZ[3], &a[3], &b[3], &c[3]);

  //something of my own to determine X,Y,Z, and rms dist
  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[1], I[1], G[1],      gX[1], gY[1], gZ[1], a[1], b[1], c[1], &x01,&y01,&z01);

  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[2], I[2], G[2],      gX[2], gY[2], gZ[2], a[2], b[2], c[2], &x02,&y02,&z02);

  point_line_line(Ex[0], I[0], G[0], mmp, gX[0], gY[0], gZ[0], a[0], b[0], c[0],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x03,&y03,&z03);

  point_line_line(Ex[1], I[1], G[1], mmp, gX[1], gY[1], gZ[1], a[1], b[1], c[1],
	              Ex[2], I[2], G[2],      gX[2], gY[2], gZ[2], a[2], b[2], c[2], &x12,&y12,&z12);

  point_line_line(Ex[1], I[1], G[1], mmp, gX[1], gY[1], gZ[1], a[1], b[1], c[1],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x13,&y13,&z13);

  point_line_line(Ex[2], I[2], G[2], mmp, gX[2], gY[2], gZ[2], a[2], b[2], c[2],
	              Ex[3], I[3], G[3],      gX[3], gY[3], gZ[3], a[3], b[3], c[3], &x23,&y23,&z23);

  *Xp=(1./6.)*(x01+x02+x03+x12+x13+x23);  
  *Yp=(1./6.)*(y01+y02+y03+y12+y13+y23);  
  *Zp=(1./6.)*(z01+z02+z03+z12+z13+z23);	
  *dist=1;

}

void det_lsq_3d (Ex, I, G, ap, mm, x1, y1, x2, y2, x3, y3, x4, y4, Xp, Yp, Zp, num_cams)
Exterior	Ex[4];
Interior	I[4];
Glass   	G[4];
ap_52		ap[4];
mm_np		mm;
double		x1, y1, x2, y2, x3, y3, x4, y4, *Xp,*Yp,*Zp;
int num_cams;
{

	    int     i,count_inner=0,n,m,flag[4];
	    double  d_inner=0.,x,y;
	    double X[4],Y[4],Z[4],a[4],b[4],c[4],dist,dist_error,X_pos[6],Y_pos[6],Z_pos[6],XX,YY,ZZ,si0,sqX,sqY,sqZ;
	    
        //new det_lsq function, bloody fast!
		flag[0]=0;flag[1]=0;flag[2]=0;flag[3]=0;
		if(x1>-999){
			flag[0]=1;
			x = x1 - I[0].xh;
	        y = y1 - I[0].yh;
	        //correct_brown_affin (x, y, ap[0], &x, &y);
		    ray_tracing_v2 (x,y, Ex[0], I[0], G[0], mmp, &X[0], &Y[0], &Z[0], &a[0], &b[0], &c[0]);
		}		
		if(x2>-999){
			flag[1]=1;
			x = x2 - I[1].xh;
	        y = y2 - I[1].yh;
	        //correct_brown_affin (x, y, ap[1], &x, &y);
		    ray_tracing_v2 (x,y, Ex[1], I[1], G[1], mmp, &X[1], &Y[1], &Z[1], &a[1], &b[1], &c[1]);
		}		
		if(x3>-999){
			flag[2]=1;
			x = x3 - I[2].xh;
	        y = y3 - I[2].yh;
	        //correct_brown_affin (x, y, ap[2], &x, &y);
		    ray_tracing_v2 (x,y, Ex[2], I[2], G[2], mmp, &X[2], &Y[2], &Z[2], &a[2], &b[2], &c[2]);
		}		
		if(x4>-999){
			flag[3]=1;
			x = x4 - I[3].xh;
	        y = y4 - I[3].yh;
	        //correct_brown_affin (x, y, ap[3], &x, &y);
		    ray_tracing_v2 (x,y, Ex[3], I[3], G[3], mmp, &X[3], &Y[3], &Z[3], &a[3], &b[3], &c[3]);
		}

		count_inner=0;
		for (n = 0; n < num_cams; n++){
			for(m = n+1; m < num_cams; m++){
				if(flag[n]==1 && flag[m]==1){
                    mid_point(X[n],Y[n],Z[n],a[n],b[n],c[n],X[m],Y[m],Z[m],a[m],b[m],c[m],&dist,&XX,&YY,&ZZ);
                    d_inner += dist;
					X_pos[count_inner]=XX;Y_pos[count_inner]=YY;Z_pos[count_inner]=ZZ;
					count_inner++;
				}
			}
		}
        d_inner/=(double)count_inner;		
		XX=0.;YY=0.;ZZ=0.;
		for(i=0;i<count_inner;i++){
           XX+=X_pos[i]; 
		   YY+=Y_pos[i];
		   ZZ+=Z_pos[i];
		}
		XX/=(double)count_inner;YY/=(double)count_inner;ZZ/=(double)count_inner;
		//end of new det_lsq
		*Xp=XX;
		*Yp=YY;
		*Zp=ZZ;

		//statistics
		si0=0.;sqX=0.;sqY=0.;sqZ=0.;
		for(i=0;i<count_inner;i++){
           si0+=pow(X_pos[i]-XX,2.)+pow(Y_pos[i]-YY,2.)+pow(Z_pos[i]-ZZ,2.);
           sqX+=pow(X_pos[i]-XX,2.); 
		   sqY+=pow(Y_pos[i]-YY,2.);
		   sqZ+=pow(Z_pos[i]-ZZ,2.);		   
		}
		si0/=(double)count_inner;sqX/=(double)count_inner;sqY/=(double)count_inner;sqZ/=(double)count_inner;
		
		mean_sigma0 += pow(si0,0.5);
        rmsX += pow(sqX,0.5);
        rmsY += pow(sqY,0.5);
        rmsZ += pow(sqZ,0.5);
		//end of statistics

}

