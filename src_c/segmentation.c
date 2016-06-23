/****************************************************************************

Routine:				segmentation.c

Author/Copyright:		Hans-Gerd Maas

Address:		       	Institute of Geodesy and Photogrammetry
				       	ETH - Hoenggerberg
				       	CH - 8093 Zurich

Creation Date:			1988/89

Description:	  target recognition with adaptive threshold
		  based on high pass filter and connectivity analysis
						with peak fitting technique

****************************************************************************/

#include "ptv.h"
#include <optv/image_processing.h>
#include "image_processing.h"


typedef struct
{
  int 	pos, status;
  short	xmin, xmax, ymin, ymax;
  int   n, sumg;
  double  x, y;
  int   unr, touch[4], n_touch;	/* unified with target unr, touching ... */
}
peak;


void targ_rec (img0, img, par_file, xmin,xmax,ymin,ymax, pix, nr, num, cpar)

unsigned char	*img, *img0;   	/* image data, image to be set to zero */
char	       	par_file[];    	/* name of parameter file */
int    		xmin,xmax,ymin,ymax;	/* search area */
target	       	pix[];		       	/* pixel coord array, global */
int	       	nr;		       	/* image number for display */
int	       	*num;	       		/* number of detections */
control_par *cpar;

/*  thresholding and center of gravity with a peak fitting technique  */
/*  uses 4 neighbours for connectivity and 8 to find local maxima     */

{
  register int	i, j, m;
  int	      	n=0, n_wait=0, n_targets=0, sumg, sumg_min;
  int	        numpix;
  int	       	thres, gvthres[4], disco;
  int	       	cr_sz;			       /* size of crosses to be drawn */
  int	       	nnmin,nnmax, nx,ny, nxmin,nxmax, nymin,nymax;
  int	       	xa,ya,xb,yb, x4[4],y4[4], xn,yn;
  double	       	x, y;
  FILE	       	*fpp;

  register unsigned char	gv, gvref;

  targpix	       	waitlist[2048];
  
  /* avoid many dereferences */
  int imx, imy;
  imx = cpar->imx;
  imy = cpar->imy;

  /* read image name, filter dimension and threshold from parameter file */
  printf("inside targ_rec (segmentation.c) \n");
  fpp = fopen(par_file,"r");
  if (fpp){
  	printf("opened file %s", &par_file);
  	printf("assigned file header %d", fpp);
  	fscanf (fpp, "%d", &gvthres[0]);      /* threshold for binarization 1.image */
  	fscanf (fpp, "%d", &gvthres[1]);      /* threshold for binarization 2.image */
  	fscanf (fpp, "%d", &gvthres[2]);      /* threshold for binarization 3.image */
  	fscanf (fpp, "%d", &gvthres[3]);      /* threshold for binarization 4.image */
  	fscanf (fpp, "%d", &disco);		/* max discontinuity */
  	fscanf (fpp, "%d  %d", &nnmin, &nnmax);	/* min. and max. number of */
  	fscanf (fpp, "%d  %d", &nxmin, &nxmax);	/* pixels per target,  	*/
  	fscanf (fpp, "%d  %d", &nymin, &nymax);	/* abs, in x, in y    	*/
  	fscanf (fpp, "%d", &sumg_min);		       	/* min. sumg */
  	fscanf (fpp, "%d", &cr_sz);				/* size of crosses */
  	fclose (fpp);
  }
  else{
  	printf("problem opening %s\n", &par_file);
  }
  /* give thres value refering to image number */
  thres=gvthres[nr];

  /*  thresholding and connectivity analysis in image  */

  for (i=ymin; i<ymax; i++)  for (j=xmin; j<xmax; j++)
    {
      gv = *(img + i*imx + j);
      if ( gv > thres)
	if (	gv >= *(img + i*imx + j-1)
		&&	gv >= *(img + i*imx + j+1)
		&&	gv >= *(img + (i-1)*imx + j)
		&&	gv >= *(img + (i+1)*imx + j)
		&&	gv >= *(img + (i-1)*imx + j-1)
		&&	gv >= *(img + (i+1)*imx + j-1)
		&&	gv >= *(img + (i-1)*imx + j+1)
		&&	gv >= *(img + (i+1)*imx + j+1) )

	  /* => local maximum, 'peak' */
	  {
	    yn=i;  xn=j;
	    xn = xn;
	    sumg = gv;  *(img + i*imx + j) = 0;
	    xa = xn;  xb = xn;  ya = yn;  yb = yn;
	    gv -= thres;
	    x = (xn) * gv;
	    y = yn * gv;
	    numpix = 1;
	    waitlist[0].x = j;  waitlist[0].y = i;  n_wait = 1;

	    while (n_wait > 0)
	      {
		gvref = *(img0 + imx*(waitlist[0].y) + (waitlist[0].x));

		x4[0] = waitlist[0].x - 1;  y4[0] = waitlist[0].y;
		x4[1] = waitlist[0].x + 1;  y4[1] = waitlist[0].y;
		x4[2] = waitlist[0].x;  y4[2] = waitlist[0].y - 1;
		x4[3] = waitlist[0].x;  y4[3] = waitlist[0].y + 1;

		for (n=0; n<4; n++)
		  {
		    xn = x4[n];  yn = y4[n];
		    xn = xn;
		    gv = *(img + imx*yn + xn);

		    /* conditions for threshold, discontinuity, image borders */
		    /* and peak fitting */
		    if (   (gv > thres)
			   && (xn>=xmin)&&(xn<xmax) && (yn>=ymin)&&(yn<ymax)
			   && (gv <= gvref+disco)
			   && (gvref + disco >= *(img0 + imx*(yn-1) + xn))
			   && (gvref + disco >= *(img0 + imx*(yn+1) + xn))
			   && (gvref + disco >= *(img0 + imx*yn + (xn-1)))
			   && (gvref + disco >= *(img0 + imx*yn + (xn+1)))  )
		      {
			sumg += gv;  *(img + imx*yn + xn) = 0;
			if (xn < xa)	xa = xn;	if (xn > xb)	xb = xn;
			if (yn < ya)	ya = yn;	if (yn > yb)	yb = yn;
			waitlist[n_wait].x = xn;	waitlist[n_wait].y = yn;
			x = x + (xn) * (gv - thres);
			y = y + yn * (gv - thres);
			numpix++;	n_wait++;
		      }
		  }

		n_wait--;
		for (m=0; m<n_wait; m++)  waitlist[m] = waitlist[m+1];
		waitlist[n_wait].x = 0;  waitlist[n_wait].y = 0;

	      }	/*  end of while-loop  */


	    /* check whether target touches image borders */
	    if (xa==xmin || ya==ymin || xb==xmax-1 || yb==ymax-1)	continue;


	    /* get targets extensions in x and y */
	    nx = xb - xa + 1;  ny = yb - ya + 1;

	    if (   (numpix >= nnmin) && (numpix <= nnmax)
		   && (nx >= nxmin) && (nx <= nxmax)
		   && (ny >= nymin) && (ny <= nymax)
		   && (sumg > sumg_min)			 )
	      {
		pix[n_targets].n = numpix;
		pix[n_targets].nx = nx;
		pix[n_targets].ny = ny;
		pix[n_targets].sumg = sumg;
		sumg -= (numpix*thres);
		x /= sumg;	x += 0.5;	y /= sumg;	y += 0.5;
		pix[n_targets].x = x;
		pix[n_targets].y = y;
		pix[n_targets].pnr = n_targets;
		pix[n_targets].tnr = -1;
        n_targets++;

		xn = x;  yn = y;
	//	drawcross (interp, (int) xn, (int) yn, cr_sz, nr, "Blue");

	      }
	  }	/*  end of if-loop  */
    }
  *num = n_targets;
}





void simple_connectivity (img0, img, par_file, xmin,xmax,ymin,ymax, pix, nr, num, cpar)
unsigned char	*img, *img0;   	/* image data, image to be set to zero */
char	       	par_file[];    	/* name of parameter file */
int	       	xmin,xmax,ymin,ymax;	/* search area */
target	       	pix[];        	/* pixel coord array, global */
int    	       	nr;    	       	/* image number for display */
int	       	*num;	       	/* number of detections */
control_par *cpar;

/*  thresholding and center of gravity with a peak fitting technique  */
/*  uses 4 neighbours for connectivity and 8 to find local maxima     */

{
  register int	i, j, m;
  int         	n=0, n_wait=0, n_targets=0, sumg, sumg_min;
  int	       	numpix;
  int	       	thres;
  int	      	cr_sz; 	      	/* size of crosses to be drawn */
  int          	nnmin,nnmax, nx,ny, nxmin,nxmax, nymin,nymax;
  int	       	xa,ya,xb,yb, x4[4],y4[4], xn,yn;
  double       	x, y;
  FILE	       	*fpp;

  register unsigned char  gv, gvref;

  targpix	waitlist[2048];

  /* avoid many dereferences */
  int imx, imy;
  imx = cpar->imx;
  imy = cpar->imy;

  /* read image name, threshold and shape limits from parameter file */
  fpp = fopen (par_file, "r");
  fscanf (fpp, "%d", &thres);    	/* threshold for binarization  	 */
  fscanf (fpp, "%d", &n);		/* threshold value for discontinuity */
  fscanf (fpp, "%d  %d", &nnmin, &nnmax);	 /* min. and max. number of  */
  fscanf (fpp, "%d  %d", &nxmin, &nxmax);	 /* pixels per target, 	*/
  fscanf (fpp, "%d  %d", &nymin, &nymax);	 /* abs, in x, in y   	*/
  fscanf (fpp, "%d", &sumg_min);	       	 /* min. sumg 		*/
  fscanf (fpp, "%d", &cr_sz);		 /* size of crosses	*/
  fclose (fpp);


  /*  thresholding and connectivity analysis in image  */

  for (i=ymin; i<ymax; i++)  for (j=xmin; j<xmax; j++)
    {
      gv = *(img + i*imx + j);
      if ( gv > 2*thres)
	{
	  yn=i;  xn=j;
	  xn = xn;
	  sumg = gv;  *(img + i*imx + j) = 0;
	  xa = xn;  xb = xn;  ya = yn;  yb = yn;
	  gv -= thres;
	  x = (xn) * gv;
	  y = yn * gv;
	  numpix = 1;
	  waitlist[0].x = j;  waitlist[0].y = i;  n_wait = 1;

	  while (n_wait > 0)
	    {
	      gvref = *(img0 + imx*(waitlist[0].y) + (waitlist[0].x));

	      x4[0] = waitlist[0].x - 1;  y4[0] = waitlist[0].y;
	      x4[1] = waitlist[0].x + 1;  y4[1] = waitlist[0].y;
	      x4[2] = waitlist[0].x;  y4[2] = waitlist[0].y - 1;
	      x4[3] = waitlist[0].x;  y4[3] = waitlist[0].y + 1;

	      for (n=0; n<4; n++)
		{
		  xn = x4[n];  yn = y4[n];
		  xn = xn;
		  gv = *(img + imx*yn + xn);

		  /* conditions for threshold, discontinuity, image borders */
		  /* and peak fitting */
		  if (   (gv > thres)
			 && (xn>=xmin)&&(xn<xmax) && (yn>=ymin)&&(yn<ymax)  )
		    {
		      sumg += gv;  *(img + imx*yn + xn) = 0;
		      if (xn < xa)	xa = xn;	if (xn > xb)	xb = xn;
		      if (yn < ya)	ya = yn;	if (yn > yb)	yb = yn;
		      waitlist[n_wait].x = xn;	waitlist[n_wait].y = yn;
		      x = x + (xn ) * (gv - thres);
		      y = y + yn * (gv - thres);
		      numpix++;	n_wait++;
		    }
		}

	      n_wait--;
	      for (m=0; m<n_wait; m++)  waitlist[m] = waitlist[m+1];
	      waitlist[n_wait].x = 0;  waitlist[n_wait].y = 0;

	    }	/*  end of while-loop  */

	  nx = xb - xa + 1;  ny = yb - ya + 1;

	  if (   (numpix >= nnmin) && (numpix <= nnmax)
		 && (nx >= nxmin) && (nx <= nxmax)
		 && (ny >= nymin) && (ny <= nymax)
		 && (sumg > sumg_min)			 )
	    {
	      pix[n_targets].n = numpix;
	      pix[n_targets].nx = nx;
	      pix[n_targets].ny = ny;
	      pix[n_targets].sumg = sumg;
	      sumg -= (numpix*thres);
	      x /= sumg;	x += 0.5;	y /= sumg;	y += 0.5;
	      pix[n_targets].x = x;
	      pix[n_targets].y = y;
	      pix[n_targets].tnr = -1;
	      pix[n_targets].pnr = n_targets;
          n_targets++;

	      xn =  x;  yn = y;
	  //    drawcross (xn, yn, cr_sz, 8);
	    }
	}	/*  end of if-loop  */
    }

  *num = n_targets;
}


