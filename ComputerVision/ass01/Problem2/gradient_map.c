#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                        Hough Transform for Circles                       */
/*                                                                          */
/*          (Copyright Andres Bruhn and Sebastian Volz, 10/2012)            */
/*                                                                          */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */

{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void dummies

     (float **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }
return;
}


/*--------------------------------------------------------------------------*/

void presmooth

     (float    **f,       /* input: original image */
      float    **u,       /* output: smoothed */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* grid size in x-direction */
      float    hy,        /* grid size in y-direction */
      float    sigma)     /* standard deviation of Gaussian */



/*
 Gaussian convolution. Copyright by Joachim Weickert 5/2000
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */


/* for sigma = 0 leave original image untouched */
 if (sigma ==0)
   {
     for (j=1; j<=ny; j++)
       for (i=1; i<=nx; i++)
   u[i][j] = f[i][j];

     return;
   }

/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma / hx) + 1;
if (length > nx)
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }


/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma / hx * sqrt(2.0 * 3.1415926))
              * exp (- (i * i * hx * hx) / (2.0 * sigma * sigma));


/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign reflecting boundary conditions */
    for (p=1; p<=length; p++)
      {
  help[length-p]      = help[length+p-1];
  help[nx+length-1+p] = help[nx+length-p];
      }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);


/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma / hy) + 1;
if (length > ny)
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma / hy * sqrt(2.0 * 3.1415927))
              * exp (- (j * j * hy * hy) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign reflecting boundary conditions */
    for (p=1; p<=length; p++)
      {
  help[length-p]      = help[length+p-1];
  help[ny+length-1+p] = help[ny+length-p];
      }

    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        u[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);


return;

} /* gauss_conv */



/*--------------------------------------------------------------------------*/

void gradient_magnitude

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x-direction */
      long    ny,          /* pixel number in y-direction */
      float   hx,          /* grid size in x-direction */
      float   hy,          /* grid size in y-direction */
      float   **u_mag)     /* gradient magnitude, output */


/*
 computes gradient magnitude
*/

{
long    i, j;       /* loop variables */
float   fx, fy;     /* local derivative in x- and y-direction */

/* mirror boundaries */
dummies (u, nx, ny);

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
        fx = (u[(long)(i+hx)][j] - u[(long)(i-hx)][j])/2;
        fy = (u[i][(long)(j+hy)] - u[i][(long)(j-hy)])/2;
        
        u_mag[i][j] = sqrt(pow(fx,2)+pow(fy,2));
     }

return;

} /* gradient_magnitude */



/*--------------------------------------------------------------------------*/

void threshold

     (float   **u,            /* image, unchanged */
      long    nx,             /* pixel number in x-direction */
      long    ny,             /* pixel number in y-direction */
      float   t_edge,         /* threshold */
      float   **u_threshold)  /* thresholded image, output */


/*
 computes thresholded image
*/

{
long    i, j;       /* loop variables */

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
       if (u[i][j]<t_edge)
   u_threshold[i][j] = 0.0;
       else
   u_threshold[i][j] = 255.0;
     }

return;

} /* threshold */


/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in[80];               /* for reading data */
char   out1[80];             /* for reading data */
float  **f;                  /* original image */
float  **u;                  /* presmoothed image */
float  **u_mag;              /* gradient magnitude */
float  **b;                  /* binary image for counting coins */
long   i, j, k;              /* loop variables */
long   nx, ny;               /* image size in x, y direction */
float  hx, hy;               /* grid size in x, y direction */
FILE   *inimage, *outimage;  /* input file, output file */
float  sigma;                /* standard deviation */
float  t_edge;               /* gradient threshold */
unsigned char byte;          /* for data conversion */


printf("\n");
printf("HOUGH TRANSFORM FOR CIRCLES\n\n");
printf("****************************************************************\n");
printf("\n");
printf("    Copyright 2012 by Andres Bruhn and Sebastian Volz   \n");
printf("    Institute for Visualization and Interactive Systems \n");
printf("    University of Stuttgart, Germany\n");
printf("\n");
printf("    All rights reserved. Unauthorized usage,\n");
printf("    copying, hiring, and selling prohibited.\n");
printf("\n");
printf("    Send bug reports to\n");
printf("    bruhn@vis.uni-stuttgart.de\n");
printf("\n");
printf("****************************************************************\n\n");

/* -------------- R E A D   I N ------------- */


/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                                            ");
gets (in);

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets(row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);

/* allocate storage */
alloc_matrix (&f, nx+2, ny+2);
alloc_matrix (&u, nx+2, ny+2);
alloc_matrix (&u_mag, nx+2, ny+2);
alloc_matrix (&b, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     f[i][j] = (float) getc (inimage);
fclose(inimage);

/* ---- read other parameters ---- */

printf("presmoothing parameter sigma:                          ");
gets(row);  sscanf(row, "%f", &sigma);
printf("gradient threshold t_edge:                             ");
gets(row);  sscanf(row, "%f", &t_edge);
printf("output image (edge image):                             ");
gets(out1);
printf("\n");


/* -------------- C O M P U T A T I O N S ------------- */

/* ---- define grid size ---- */
hx = 1.0;
hy = 1.0;

/* ---- presmooth image ---- */
presmooth (f, u, nx, ny, hx, hy, sigma);

/* ---- compute gradient magnitude ---- */
gradient_magnitude (u, nx, ny, hy, hy, u_mag);

/* ---- threshold gradient magnitude ---- */
threshold (u_mag, nx, ny, t_edge, u_mag);

/* -------------- W R I T E   O U T ------------- */

/* ---- write output image (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out1, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# Edge Map \n");
fprintf (outimage, "# initial image:         %s\n", in);
fprintf (outimage, "# presmoothing sigma:    %f\n", sigma);
fprintf (outimage, "# edge threshold:        %f\n", t_edge);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     if (u_mag[i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (u_mag[i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(u_mag[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out1);

/* ---- disallocate storage ---- */
disalloc_matrix (f, nx+2, ny+2);
disalloc_matrix (u, nx+2, ny+2);
disalloc_matrix (u_mag, nx+2, ny+2);
disalloc_matrix (b, nx+2, ny+2);

return(0);
}
