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

void unsharp_mask (n, img0, img_lp)

int          n;
unsigned char *img0, *img_lp;


{
   register unsigned char  *imgum, *ptrl, *ptrr, *ptrz;
   int                 *buf1, *buf2, buf, *end;
   register int            *ptr, *ptr1, *ptr2, *ptr3;
   int                 ii, n2, nq, m;
   register int            i;
   //imgsize=len1; // denis
   n2 = 2*n + 1;  nq = n2 * n2;
   printf("inside unsharp_mask\n");

   imgum = (unsigned char *) calloc (imgsize, 1);
   if ( ! imgum)
   {
       puts ("calloc for imgum --> error");  exit (1);
   }
   printf("after calloc unsharp_mask\n");      
   buf1 = (int *) calloc (imgsize, sizeof(int));
   if ( ! buf1)
   {
       puts ("calloc for buf1 --> error");  exit (1);
   }
   buf2 = (int *) calloc (cpar->imx, sizeof(int));

   printf("after calloc2 unsharp_mask\n"); 

   /* set imgum = img0 (so there cannot be written to the original image) */
   for (ptrl=imgum, ptrr=img0; ptrl<(imgum+imgsize); ptrl++, ptrr++)
   {
     *ptrl = *ptrr;


   }   

   /* cut off high gray values (not in general use !)
   for (ptrz=imgum; ptrz<(imgum+imgsize); ptrz++) if (*ptrz > 160) *ptrz = 160; */




   /* --------------  average over lines  --------------- */

   for (i=0; i < cpar->imy; i++)
   {
       ii = i * cpar->imx;
       /* first element */
       buf = *(imgum+ii);  *(buf1+ii) = buf * n2;
       
       /* elements 1 ... n */
       for (ptr=buf1+ii+1, ptrr=imgum+ii+2, ptrl=ptrr-1, m=3;
            ptr<buf1+ii+n+1; ptr++, ptrl+=2, ptrr+=2, m+=2)
       {
           buf += (*ptrl + *ptrr);
           *ptr = buf * n2 / m;
       }
       
       /* elements n+1 ... imx-n */
       for (ptrl=imgum+ii, ptr=buf1+ii+n+1, ptrr=imgum+ii+n2;
            ptrr<imgum+ii+cpar->imx; ptrl++, ptr++, ptrr++)
       {
           buf += (*ptrr - *ptrl);
           *ptr = buf;
       }
       
       /* elements imx-n ... imx */
       for (ptrl=imgum+ii+cpar->imx-n2, ptrr=ptrl+1, ptr=buf1+ii+cpar->imx-n, m=n2-2;
            ptr<buf1+ii+cpar->imx; ptrl+=2, ptrr+=2, ptr++, m-=2)
       {
           buf -= (*ptrl + *ptrr);
           *ptr = buf * n2 / m;
       }
   }
   
   free (imgum);


   /* -------------  average over columns  -------------- */

   end = buf2 + cpar->imx;

   /* first line */
   for (ptr1=buf1, ptr2=buf2, ptrz=img_lp; ptr2<end; ptr1++, ptr2++, ptrz++)
   {
       *ptr2 = *ptr1;
       *ptrz = *ptr2/n2;
   }
   
   /* lines 1 ... n */
   for (i=1; i<n+1; i++)
   {
       ptr1 = buf1 + (2*i-1)*cpar->imx;
       ptr2 = ptr1 + cpar->imx;
       ptrz = img_lp + i*cpar->imx;
       for (ptr3=buf2; ptr3<end; ptr1++, ptr2++, ptr3++, ptrz++)
       {
           *ptr3 += (*ptr1 + *ptr2);
           *ptrz = n2 * (*ptr3) / nq / (2*i+1);
       }
   }
   
   /* lines n+1 ... imy-n-1 */
   for (i=n+1, ptr1=buf1, ptrz=img_lp+cpar->imx*(n+1), ptr2=buf1+cpar->imx*n2;
        i<cpar->imy-n; i++)
   {
       for (ptr3=buf2; ptr3<end; ptr3++, ptr1++, ptrz++, ptr2++)
       {
           *ptr3 += (*ptr2 - *ptr1);
           *ptrz = *ptr3/nq;
       }
   }
   
   /* lines imy-n ... imy */
   for (i=n; i>0; i--)
   {
       ptr1 = buf1 + (cpar->imy-2*i-1)*cpar->imx;
       ptr2 = ptr1 + cpar->imx;
       ptrz = img_lp + (cpar->imy-i)*cpar->imx;
       for (ptr3=buf2; ptr3<end; ptr1++, ptr2++, ptr3++, ptrz++)
       {
           *ptr3 -= (*ptr1 + *ptr2);
           *ptrz = n2 * (*ptr3) / nq / (2*i+1);
       }
   }
   
   
   free (buf1);
   printf("end unsharp_mask\n");  
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
