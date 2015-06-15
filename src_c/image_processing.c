/****************************************************************************

Routine:	       	image_processing.c

Author/Copyright:      	Hans-Gerd Maas

Address:	       	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	1988
	
Description:	       	different image processing routines ...
	
Routines contained:    	filter_3:	3*3 filter, reads matrix from filter.par
		       	lowpass_3:	3*3 local average with 9 pointers, fast
		       	lowpass_n:	n*n local average,, fast
		       			computation time independent from n
		       	histogram:	computes histogram
		       	enhance:	enhances gray value spectrum to 0..255,
		      			some extreme gray values are cut off
		       	mark_img:	reads image and pixel coordinate set,
			       		marks the image in a certain (monocomp)
			       		color and writes image

****************************************************************************/

#include "ptv.h"


void filter_3(unsigned char *img, unsigned char *img_lp, control_par *cpar) {
	register unsigned char	*ptr, *ptr1, *ptr2, *ptr3,
		             	*ptr4, *ptr5, *ptr6,
	                        *ptr7, *ptr8, *ptr9;
	int	       	end;
	float	       	m[9], sum;
	short	       	buf;
	register int	i;
	FILE	       	*fp;

	/* read filter elements from parameter file */
	fp = fopen_r ("filter.par");
	for (i=0, sum=0; i<9; i++)
	{
		fscanf (fp, "%f", &m[i]);
		sum += m[i];
	}
	fclose (fp);  if (sum == 0) exit(1);
	
	end = imgsize - 513;
	
	ptr  = img_lp + 513;
	ptr1 = img;
    ptr2 = img + 1;
    ptr3 = img + 2;
    
	ptr4 = img + cpar->imx;
    ptr5 = ptr4 + 1;
    ptr6 = ptr4 + 2;
    
	ptr7 = img + 2*cpar->imx;
    ptr8 = ptr7 + 1;
    ptr9 = ptr7 + 2;

	for (i=513; i<end; i++)
	{
		buf = m[0] * *ptr1++  +  m[1] * *ptr2++  +  m[2] * *ptr3++
			+ m[3] * *ptr4++  +  m[4] * *ptr5++  +  m[5] * *ptr6++
			+ m[6] * *ptr7++  +  m[7] * *ptr8++  +  m[8] * *ptr9++;
		buf /= sum;    if (buf > 255)  buf = 255;    if (buf < 8)  buf = 8;
		*ptr++ = buf;
	}
}

void histogram (img, hist)

unsigned char	*img;
int	       	*hist;

{
	int	       	i;
	unsigned char  	*end;
	register unsigned char	*ptr;

	
	for (i=0; i<256; i++)  hist[i] = 0;
	
	end = img + imgsize;
	for (ptr=img; ptr<end; ptr++)
	{
		hist[*ptr]++;
	}
}






void lowpass_3(unsigned char *img, unsigned char *img_lp, control_par *cpar) {
	register unsigned char	*ptr,*ptr1,*ptr2,*ptr3,*ptr4,
		       		*ptr5,*ptr6,*ptr7,*ptr8,*ptr9;
	short  	       		buf;
	register int   		i;
	
	ptr  = img_lp + 513;
	ptr1 = img;
	ptr2 = img + 1;
	ptr3 = img + 2;
    
	ptr4 = img + cpar->imx;
	ptr5 = ptr4 + 1;
	ptr6 = ptr4 + 2;
    
	ptr7 = img + 2*cpar->imx;
	ptr8 = ptr7 + 1;
	ptr9 = ptr7 + 2;

	for (i=0; i<imgsize; i++)
	{
		buf = *ptr5++ + *ptr1++ + *ptr2++ + *ptr3++ + *ptr4++
					  + *ptr6++ + *ptr7++ + *ptr8++ + *ptr9++;
		*ptr++ = buf/9;
	}

}


void split (img, field, cpar)

unsigned char	*img;
int	       	field;
control_par *cpar;

{
	register int   		i, j;
	register unsigned char	*ptr;
	unsigned char	       	*end;

	switch (field)
	{
		case 0:  /* frames */
				return;	 break;

        case 1:  /* odd lines */
            for (i = 0; i < cpar->imy/2; i++)  
                for (j = 0; j < cpar->imx; j++)
                    *(img + cpar->imx*i + j) = *(img + 2*cpar->imx*i + j + cpar->imx);  
                    break;

        case 2:  /* even lines */
            for (i = 0; i < cpar->imy/2; i++)
                for (j = 0; j < cpar->imx; j++)
                    *(img + cpar->imx*i + j) = *(img + 2*cpar->imx*i + j);  break;
	}
	
	end = img + imgsize;
	for (ptr=img+imgsize/2; ptr<end; ptr++)  *ptr = 2;
}






void copy_images (img1, img2)

unsigned char	*img1, *img2;

{
	register unsigned char 	*ptr1, *ptr2;
	unsigned char	       	*end;


	for (end=img1+imgsize, ptr1=img1, ptr2=img2; ptr1<end; ptr1++, ptr2++)
	*ptr2 = *ptr1;
}



/*------------------------------------------------------------------------
	Subtract mask, Matthias Oswald, Juli 08
  ------------------------------------------------------------------------*/
void subtract_mask (img, img_mask, img_new) 

unsigned char	*img, *img_mask, *img_new;

{
	register unsigned char 	*ptr1, *ptr2, *ptr3;
	int i;
	
	for (i=0, ptr1=img, ptr2=img_mask, ptr3=img_new; i<imgsize; ptr1++, ptr2++, ptr3++, i++)
    {
      if (*ptr2 == 0)  *ptr3 = 0;
      else  *ptr3 = *ptr1;
    }
 }
