#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                   OPTIC FLOW USING LUCAS AND KANADE                      */
/*                                                                          */
/*                 (Copyrights Joachim Weickert, 2/2003,                    */
/*                  and Stephan Didas, 1/2007)                              */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 Optic flow calculation with method by Lucas and Kanade
*/


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

/* ---------------------------------------------------------------------- */

void gauss_conv

     (float    sigma,     /* standard deviation of Gaussian */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    precision, /* cutoff at precision * sigma */
      long     bc,        /* type of boundary condition */
                          /* 0=Dirichlet, 1=reflecing, 2=periodic */
      float    **f)       /* input: original image ;  output: smoothed */


/*
 Gaussian convolution.
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */


/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(precision * sigma / hx) + 1;
if ((bc != 0) && (length > nx))
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926))
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

    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[nx+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[nx+length-1+p] = help[nx+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[nx+length-p];
           help[nx+length-1+p] = help[length+p-1];
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
length = (long)(precision * sigma / hy) + 1;
if ((bc != 0) && (length > ny))
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927))
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
    /* assign boundary conditions */
    if (bc == 0) /* Dirichlet boundary conditions */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = 0.0;
           help[ny+length-1+p] = 0.0;
           }
    else if (bc == 1) /* reflecting b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[length+p-1];
           help[ny+length-1+p] = help[ny+length-p];
           }
    else if (bc == 2) /* periodic b.c. */
       for (p=1; p<=length; p++)
           {
           help[length-p]      = help[ny+length-p];
           help[ny+length-1+p] = help[length+p-1];
           }

    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        f[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);

return;

} /* gauss_conv */

/* ---------------------------------------------------------------------- */

void create_eq_systems

     (float    **f1,      /* frame 1, input */
      float    **f2,      /* frame 2, input */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    hx,        /* pixel size in x direction */
      float    hy,        /* pixel size in y direction */
      float    ht,        /* distance between two frames */
      float    rho,       /* integration scale */
      float    **dxx,     /* element of structure tensor, output */
      float    **dxy,     /* element of structure tensor, output */
      float    **dyy,     /* element of structure tensor, output */
      float    **dxz,     /* element of the right-hand side, output */
      float    **dyz)     /* element of the right-hand side, output */

/*
 Calculates the entries of the structure tensor and the right-hand side
 for the linear systems of equations.
*/

{
long    i, j;                 /* loop variables */
float   df_dx, df_dy, df_dz;  /* derivatives of f */
float   w1, w2, w3, w4, w5;   /* time savers */


/* ---- calculate gradient and its tensor product ---- */

dummies(f1, nx, ny); 
dummies(f2, nx, ny); 

for(i=1; i<=nx; i++) {
    for(j=1; j<= ny; j++) {
	 /* compute the spatial derivatives using Sobel operators */
	 w1 = 1.0 / (8.0 * hx);
	 w2 = 1.0 / (4.0 * hx);
	 
	 df_dx = w1 * (f1[i+1][j+1] - f1[i-1][j+1]
		       + f1[i+1][j-1] - f1[i-1][j-1])
	       + w2 * (f1[i+1][j] - f1[i-1][j]);
	 
	 w3 = 1.0 / (8.0 * hy);
	 w4 = 1.0 / (4.0 * hy);
	 
	 df_dy = w3 * (f1[i+1][j+1] - f1[i+1][j-1]
		       + f1[i-1][j+1] - f1[i-1][j-1])
               + w4 * (f1[i][j+1] - f1[i][j-1]);

	 /* time derivative with forward difference, ht = 1 */
	 w5 = 1.0 / ht;
	 df_dz = w5 * (f2[i][j] - f1[i][j]);
	 
	 /* calculate matrix entries and right-hand side */
	 dxx[i][j] = df_dx * df_dx;
	 dxy[i][j] = df_dx * df_dy;
	 dyy[i][j] = df_dy * df_dy;
	 dxz[i][j] = -1 * df_dx * df_dz;
	 dyz[i][j] = -1 * df_dy * df_dz;
     }
 }
/* ---- smoothing at integration scale, Dirichlet b.c. ---- */

if (rho > 0.0)
   {
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxx);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxy);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dyy);
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dxz):
   gauss_conv (rho, nx, ny, hx, hy, 5.0, 0, dyz);
   }

return;
} /* create_eq_systems */

/* ---------------------------------------------------------------------- */

void lucas_kanade 

     (long     nx,          /* image dimension in x direction */ 
      long     ny,          /* image dimension in y direction */ 
      float    **dxx,       /* component of system matrix, input */
      float    **dxy,       /* component of system matrix, input */
      float    **dyy,       /* component of system matrix, input */
      float    **dxz,       /* component of right-hand side, input */
      float    **dyz,       /* component of right-hand side, input */
      float    eps,         /* threshold */
      float    **c,         /* confidence of the flow, output */
      float    **u,         /* x component of optic flow, output */
      float    **v)         /* v component of optic flow, output */


/* 
   Performs optic flow estimation with Lucas-Kanade method.
*/

{
long i, j;       /* loop variables */
float trace;     /* trace of the matrix at pixel [i][j] */
float det;       /* determinant of the matrix at pixel [i][j] */

for(i=1; i<=nx; i++) {
    for(j=1; j<=ny; j++) {

	u[i][j] = (dxy[i][j]*dyz[i][j] - dxz[i][j]*dyy[i][j])/
		(dxx[i][j]*dyy[i][j] - dxy[i][j]*dxy[i][j]);
	
	u[i][j] = (dxz[i][j]*dxy[i][j] - dyz[i][j]*dxx[i][j])/
		(dxx[i][j]*dyy[i][j] - dxy[i][j]*dxy[i][j]);
    }
}

return;
} /* lucas_kanade */

/* ---------------------------------------------------------------------- */

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *vari)       /* variance, output */

/*
 calculates minimum, maximum, mean and variance of an image u
*/

{
long    i, j;       /* loop variables */
float   help;       /* auxiliary variable */
double  help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help2 = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help2 = help2 + (double)u[i][j];
     }
*mean = (float)help2 / (nx * ny);

*vari = 0.0;
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     {
     help  = u[i][j] - *mean;
     *vari = *vari + help * help;
     }
*vari = *vari / (nx * ny);

return;

} /* analyse */

/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in1[80], in2[80];     /* for reading data */
char   out1[80], out2[80];   /* for reading data */
float  **f1, **f2;           /* images */
float  **dxx, **dxy, **dyy;  /* components of matrix */
float  **dxz, **dyz;         /* components of right-hand side */
float  **u, **v;             /* optic flow components */
float  **w, **c;             /* optic flow magnitude and confidence */
long   i, j;                 /* loop variables */
long   nx, ny;               /* image size in x, y direction */
FILE   *inimage, *outimage;  /* input file, output file */
float  rho;                  /* integration scale */
float  eps;                  /* threshold for Lucas-Kanade method */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  vari;                 /* variance */
float  help;                 /* time saver */
unsigned char byte;          /* for data conversion */


printf("\n");
printf("OPTIC FLOW COMPUTATION WITH THE METHOD OF LUCAS AND KANADE\n\n");
printf("***********************************************\n\n");
printf("  Copyright by Joachim Weickert 2/2003 and     \n");
printf("  Stephan Didas 1/2007                         \n");
printf("  Faculty of Mathematics and Computer Science  \n");
printf("  Saarland University, Saarbruecken, Germany.  \n\n");
printf("  All rights reserved. Unauthorized usage,     \n");
printf("  copying, hiring, and selling prohibited.     \n\n");
printf("  Send bug reports to                          \n");
printf("  bruhn@vis.uni-stuttgart.de                   \n");
printf("***********************************************\n\n");


/* ---- read input image f1 (pgm format P5) ---- */

/* read image name */
printf("input image 1:                        ");
gets (in1);

/* open pgm file and read header */
inimage = fopen(in1,"r");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets(row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);

/* allocate storage */
alloc_matrix (&f1, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     f1[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read input image f2 (pgm format P5) ---- */

/* read image name */
printf("input image 2:                        ");
gets (in2);

/* open pgm file and read header */
inimage = fopen(in2,"r");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets(row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);

/* allocate storage */
alloc_matrix (&f2, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     f2[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read other parameters ---- */

printf("integration scale rho (>0):           ");
gets(row);  sscanf(row, "%f", &rho);
printf("threshold epsilon (>0):               ");
gets(row);  sscanf(row, "%f", &eps);
printf("output image (flow classification):   ");
gets(out1);
printf("output image (flow magnitude):        ");
gets(out2);
printf("\n");


/* ---- initializations ---- */
/* allocate storage */
alloc_matrix (&dxx, nx+2, ny+2);
alloc_matrix (&dxy, nx+2, ny+2);
alloc_matrix (&dyy, nx+2, ny+2);
alloc_matrix (&dxz, nx+2, ny+2);
alloc_matrix (&dyz, nx+2, ny+2);
alloc_matrix (&u, nx+2, ny+2);
alloc_matrix (&v, nx+2, ny+2);
alloc_matrix (&w, nx+2, ny+2);
alloc_matrix (&c, nx+2, ny+2);

create_eq_systems(f1, f2, nx, ny, 1.0, 1.0, 1.0, rho, dxx, dxy, dyy, 
                  dxz, dyz);

lucas_kanade(nx, ny, dxx, dxy, dyy, dxz, dyz, eps, c, u, v);

/* calculate flow magnitude */
for (i=1; i<=nx; i++)
    for (j=1; j<=ny; j++)
	w[i][j] = sqrt (u[i][j] * u[i][j] + v[i][j] * v[i][j]);

/* check minimum, maximum, mean, variance of flow magnitude */
analyse (w, nx, ny, &min, &max, &mean, &vari);
printf("maximum:       %8.2f \n", max);
printf("mean:          %8.2f \n", mean);
printf("variance:      %8.2f \n", vari);

/* ---- write output image c (pgm format P5) ---- */

outimage = fopen (out1, "w");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# optic flow confidence, Lucas and Kanade method\n");
fprintf (outimage, "# initial image 1:  %s\n", in1);
fprintf (outimage, "# initial image 2:  %s\n", in2);
fprintf (outimage, "# rho:              %8.4f\n", rho);
fprintf (outimage, "# threshold eps:    %8.4f\n", eps);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     help = c[i][j];
     if (help < 0.0)
        byte = (unsigned char)(0.0);
     else if (help > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(help);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);

printf("\n");
printf("output image %s successfully written\n\n", out1);


/* ---- write output image w (pgm format P5) ---- */

outimage = fopen (out2, "w");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# optic flow magnitude, Lucas and Kanade method\n");
fprintf (outimage, "# initial image 1:  %s\n", in1);
fprintf (outimage, "# initial image 2:  %s\n", in2);
fprintf (outimage, "# rho:              %8.4f\n", rho);
fprintf (outimage, "# threshold eps:    %8.4f\n", eps);
fprintf (outimage, "# max:              %8.4f\n", max);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     help = 255.0 * w[i][j] / max;
     if (help < 0.0)
        byte = (unsigned char)(0.0);
     else if (help > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(help);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);

printf("\n");
printf("output image %s successfully written\n\n", out2);


/* ---- disallocate storage ---- */
disalloc_matrix (f1, nx+2, ny+2);
disalloc_matrix (f2, nx+2, ny+2);
disalloc_matrix (u,  nx+2, ny+2);
disalloc_matrix (v,  nx+2, ny+2);
disalloc_matrix (w,  nx+2, ny+2);
disalloc_matrix (dxx, nx+2, ny+2);
disalloc_matrix (dxy, nx+2, ny+2);
disalloc_matrix (dyy,  nx+2, ny+2);
disalloc_matrix (dxz,  nx+2, ny+2);
disalloc_matrix (dyz,  nx+2, ny+2);
disalloc_matrix (c,  nx+2, ny+2);

return(0);
}
