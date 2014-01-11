#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libraries.h"

/*************************/
/* PDE morph
 */


/*--------------------------------------------------------------------------*/
inline float max(float a, float b) {return (a>b)?a:b;}
/*--------------------------------------------------------------------------*/
inline float min(float a, float b) {return (a<b)?a:b;}
/*--------------------------------------------------------------------------*/
inline float sq(float a) {return a*a;}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

inline float dilation_point
  (
  float **u,
  float tau,
  float hx_1,
  float hy_1,
  int i,
  int j
  )
  {
  return ((sqrt( (
	sq(min( (u[i][j] - u[i-1][j])/hx_1, 0)) + 
	sq(max( (u[i+1][j] - u[i][j])/hx_1, 0)) + 
	sq(min( (u[i][j] - u[i][j-1])/hy_1, 0)) +
	sq(max( (u[i][j+1] - u[i][j])/hy_1, 0))) )
	* tau) + u[i][j]);
  }

/*--------------------------------------------------------------------------*/

/* dilation iteration step */
void dilation
  (
  float **u,
  int nx,
  int ny,
  float hx,
  float hy,
  float tau
  )
{
int i;                        // loop variables
int j;                        // loop variables

float hx_1;                   // time savers
float hy_1;                   // time savers

float **tmp;                  // temporary array

/* compute time savers */
hx_1 = 1.0 / hx;
hy_1 = 1.0 / hy;

/* allocate temporary memory */
alloc_matrix(&tmp,nx+2,ny+2);

/* iteration step */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  tmp[i][j] = dilation_point(u,tau,hx_1,hy_1,i,j);
  }

/* store data */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  u[i][j] = tmp[i][j];

/* free temporary memory */
disalloc_matrix(tmp,nx+2,ny+2);

return;
}

/*--------------------------------------------------------------------------*/

inline float erosion_point
  (
  float **u,
  float tau,
  float hx_1,
  float hy_1,
  int i,
  int j
  )
  {
  return ((-sqrt( (
	sq(max( (u[i][j] - u[i-1][j])/hx_1, 0)) + 
	sq(min( (u[i+1][j] - u[i][j])/hx_1, 0)) + 
	sq(max( (u[i][j] - u[i][j-1])/hy_1, 0)) +
	sq(min( (u[i][j+1] - u[i][j])/hy_1, 0))) )
	* tau) + u[i][j]);
  }

/*--------------------------------------------------------------------------*/

/* erosion iteration step */
void erosion
  (
  float **u,
  int nx,
  int ny,
  float hx,
  float hy,
  float tau
  )
{
int i;                        // loop variables
int j;                        // loop variables

float hx_1;                   // time savers
float hy_1;                   // time savers

float **tmp;                  // temporary array

/* compute time savers */
hx_1 = 1.0 / hx;
hy_1 = 1.0 / hy;

/* allocate temporary memory */
alloc_matrix(&tmp,nx+2,ny+2);

/* iteration step */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  tmp[i][j] = erosion_point(u,tau,hx_1,hy_1,i,j);
  }

/* store data */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  u[i][j] = tmp[i][j];

/* free temporary memory */
disalloc_matrix(tmp,nx+2,ny+2);

return;
}

/*--------------------------------------------------------------------------*/

/* compute structure tensor (joint strucure tensor if color image) */
void structure_tensor
  (
  float **u1,                 // input image channel 1
  float **u2,                 // input image channel 2
  float **u3,                 // input image channel 3
  float **J11,                // structure tensor entry 11
  float **J12,                // structure tensor entry 12
  float **J22,                // structure tensor entry 22
  int nx,                     // image size in x-direction
  int ny,                     // image size in y-direction
  float hx,                   // grid spacing in x-direction
  float hy,                   // grid spacing in y-direction
  float rho,                  // structure tensor convolution scale
  int color                   // switch for color / gray images
  )
{
int i;                        // loop variables
int j;                        // loop variables

float hx_1;                   // time savers
float hy_1;                   // time savers

float ux;                     // time savers
float uy;                     // time savers

/* mirror boundaries */
dummies(u1,nx,ny);
/* u2 and u3 for color images only */
if (color) dummies(u2,nx,ny);
if (color) dummies(u3,nx,ny);

/* compute time savers*/
hx_1 = 1.0 / ( 2.0 * hx );
hy_1 = 1.0 / ( 2.0 * hy );

for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  /* compute derivatives */
  ux = ( u1[i+1][j] - u1[i-1][j] ) * hx_1;
  uy = ( u1[i][j+1] - u1[i][j-1] ) * hy_1;

  /* set up structure tensor */
  J11[i][j] = ux * ux;
  J12[i][j] = ux * uy;
  J22[i][j] = uy * uy;

  /* for color images, repeat for other channels and sum up */
  if (color)
    {
    ux = ( u2[i+1][j] - u2[i-1][j] ) * hx_1;
    uy = ( u2[i][j+1] - u2[i][j-1] ) * hy_1;

    J11[i][j] += ux * ux;
    J12[i][j] += ux * uy;
    J22[i][j] += uy * uy;

    ux = ( u3[i+1][j] - u3[i-1][j] ) * hx_1;
    uy = ( u3[i][j+1] - u3[i][j-1] ) * hy_1;

    J11[i][j] += ux * ux;
    J12[i][j] += ux * uy;
    J22[i][j] += uy * uy;
    }
  }


if (rho>0.0)
  {
  /* convolve structure tensor entries */
  presmooth(J11,J11,nx,ny,rho);
  presmooth(J12,J12,nx,ny,rho);
  presmooth(J22,J22,nx,ny,rho);
  }

return;
} 

/*--------------------------------------------------------------------------*/

/* compute the eigenvector of J that belongs to the largest eigenvalue */
void get_directions

     (float **J11,        /* coeffs of (2*2)-matrix */
      float **J12,        /* coeffs of (2*2)-matrix */
      float **J22,        /* coeffs of (2*2)-matrix */
      int nx,
      int ny,
      float **r1,        /* 1. comp. of 1. eigenvector, output */
      float **r2)        /* 2. comp. of 1. eigenvector, output */



{
int i,j;
float  help, norm;    /* time savers */

for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  /* eigenvalues */
  help  = sqrt (powf (J22[i][j]-J11[i][j], 2.0) + 4.0 * J12[i][j] * J12[i][j]);

  if (help == 0.0)
     /* isotropic situation, eigenvectors arbitrary */
     {
     r1[i][j] = 1.0;
     r2[i][j] = 0.0;
     }

  else if (J11[i][j] > J22[i][j])
     {
     r1[i][j] = J11[i][j] - J22[i][j] + help;
     r2[i][j] = 2.0 * J12[i][j];
     }

  else
     {
     r1[i][j] = 2.0 * J12[i][j];
     r2[i][j] = J22[i][j] - J11[i][j] + help;
     }


  /* ---- normalize eigenvectors ---- */

  norm = sqrt (r1[i][j] * r1[i][j] + r2[i][j] * r2[i][j]);
  if (norm > 0.0)
     {
     r1[i][j] = r1[i][j] / norm;
     r2[i][j] = r2[i][j] / norm;
     }
  else
     {
     r1[i][j] = 1.0;
     r2[i][j] = 0.0;
     }  
  }

return;

} /* get_directions */
    
/*--------------------------------------------------------------------------*/

/* compute hessian matrix (joint hessian matrix if color image) */
void hessian
  (
  float **v1,                 // presmoothed image, channel 1
  float **v2,                 // presmoothed image, channel 2
  float **v3,                 // presmoothed image, channel 3
  float **H11,                // hessian matrix, entry 11
  float **H12,                // hessian matrix, entry 12
  float **H22,                // hessian matrix, entry 22
  int nx,                     // image size in x-direction
  int ny,                     // image size in y-direction
  float hx,                   // grid spacing in x-direction
  float hy,                   // grid spacing in y-direction
  int color                   // switch for color / gray images
  )
{
int i;                        // loop variables
int j;                        // loop variables

float hx_2;                   // time savers
float hy_2;                   // time savers
float hxy;                    // time savers

/* compute time savers */
hx_2 = 1.0 / ( hx * hx );
hy_2 = 1.0 / ( hy * hy );

hxy  = 1.0 / ( 4.0 * hx * hy );

for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  /*! TODO !*/
  /*! Supplement missing code here !*/
  /*! i.e. compute the hessian matrix !*/
  /*! for handling color images, check the routine structure_tensor !*/
  /*! !*/
  }

return;
}

/*--------------------------------------------------------------------------*/

/* shock filter iteration step */
void shock_filter
  (
  float **u1,                 // input image channel 1
  float **u2,                 // input image channel 2
  float **u3,                 // input image channel 3
  int nx,                     // image size in x-direction
  int ny,                     // image size in y-direction
  float hx,                   // grid spacing in x-direction
  float hy,                   // grid spacing in y-direction
  float tau,                  // time step for iterative scheme
  float sigma,                // presmoothing scale (hessian matrix)
  float rho,                  // convolution scale (structure tensor)
  int color                   // switch for color / gray images
  )
{
int i;                        // loop variables
int j;                        // loop variables

float hx_1;                   // time savers
float hy_1;                   // time savers

float **J11;                  // structure tensor, entry 11
float **J12;                  // structure tensor, entry 12
float **J22;                  // structure tensor, entry 22

float **r1;                   // 1st eigenvector of structure tensor, 1st comp.
float **r2;                   // 1st eigenvector of structure tensor, 2nd comp.

float **v1;                   // image, presmoothed (hessian)
float **v2;                   // image, presmoothed (hessian)
float **v3;                   // image, presmoothed (hessian)

float **H11;                  // hessian matrix, entry 11
float **H12;                  // hessian matrix, entry 12
float **H22;                  // hessian matrix, entry 22

float **tmp1;                 // temporary array
float **tmp2;                 // temporary array
float **tmp3;                 // temporary array

float v_eta_eta;              // second derivative in direction (r_1,r_2)^T

/* allocate memory */
alloc_matrix(&J11,nx+2,ny+2);
alloc_matrix(&J12,nx+2,ny+2);
alloc_matrix(&J22,nx+2,ny+2);

alloc_matrix(&H11,nx+2,ny+2);
alloc_matrix(&H12,nx+2,ny+2);
alloc_matrix(&H22,nx+2,ny+2);

alloc_matrix(&tmp1,nx+2,ny+2);
/* tmp2 and tmp3 for color images only */
if (color) alloc_matrix(&tmp2,nx+2,ny+2);
if (color) alloc_matrix(&tmp3,nx+2,ny+2);

alloc_matrix(&r1, nx+2,ny+2);
alloc_matrix(&r2, nx+2,ny+2);

alloc_matrix(&v1 ,nx+2,ny+2);
/* tmp2 and tmp3 for color images only */
if (color) alloc_matrix(&v2 ,nx+2,ny+2);
if (color) alloc_matrix(&v3 ,nx+2,ny+2);

/* compute structure tensor */
structure_tensor(u1,u2,u3,J11,J12,J22,nx,ny,hx,hy,rho,color);

/* get_dominant directions */
get_directions(J11,J12,J22,nx,ny,r1,r2);

/* presmooth image */
presmooth(u1,v1,nx,ny,sigma);
/* u2 and u3 for color images only */
if (color) presmooth(u2,v2,nx,ny,sigma);
if (color) presmooth(u3,v3,nx,ny,sigma);

/* mirror boundaries */
dummies(v1,nx,ny);
/* v2 and v3 for color images only */
if (color) dummies(v2,nx,ny);
if (color) dummies(v3,nx,ny);

/* compute hessian matrix */
hessian(v1,v2,v3,H11,H12,H22,nx,ny,hx,hy,color);

/* compute time_savers */
hx_1 = 1.0 / hx;
hy_1 = 1.0 / hy;



/* iteration step */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  /* compute directional derivative */
  /*! TODO !*/
  /*! Supplement missing code here !*/
  /*! i.e. compute the second order directional derivative !*/
  /*! !*/

  /* update this pixel */
  /*! TODO !*/
  /*! Supplement missing code here !*/
  /*! i.e. assign "v_eta_eta > 0.0"and "v_eta_eta < 0.0" correctly !*/
  /*! !*/
  if ( 0/*! TODO !*/ )
    {
    /* dilation: upwind scheme */
    tmp1[i][j] = dilation_point(u1,tau,hx_1,hy_1,i,j);
    /* u2 and u3 for color images only */
    if (color) tmp2[i][j] = dilation_point(u2,tau,hx_1,hy_1,i,j);
    if (color) tmp3[i][j] = dilation_point(u3,tau,hx_1,hy_1,i,j);
    }
  else if ( 1/*! TODO !*/ )
    {
    /* erosion: upwind scheme */
    tmp1[i][j] = erosion_point(u1,tau,hx_1,hy_1,i,j);
    /* u2 and u3 for color images only */
    if (color) tmp2[i][j] = erosion_point(u2,tau,hx_1,hy_1,i,j);
    if (color) tmp3[i][j] = erosion_point(u3,tau,hx_1,hy_1,i,j);
    }
  else if ( v_eta_eta == 0.0 )
    {
    /* none of them */
    tmp1[i][j] = u1[i][j];
    /* u2 and u3 for color images only */
    if (color) tmp2[i][j] = u2[i][j];
    if (color) tmp3[i][j] = u3[i][j];
    }
  }


  
/* copy data */
for (i=1;i<=nx;i++)
for (j=1;j<=ny;j++)
  {
  u1[i][j] = tmp1[i][j];
    /* u2 and u3 for color images only */
  if (color) u2[i][j] = tmp2[i][j];
  if (color) u3[i][j] = tmp3[i][j];
  }

/* free memory */
disalloc_matrix(J11,nx+2,ny+2);
disalloc_matrix(J12,nx+2,ny+2);
disalloc_matrix(J22,nx+2,ny+2);

disalloc_matrix(H11,nx+2,ny+2);
disalloc_matrix(H12,nx+2,ny+2);
disalloc_matrix(H22,nx+2,ny+2);

disalloc_matrix(tmp1,nx+2,ny+2);
/* tmp2 and tmp3 for color images only */
if (color) disalloc_matrix(tmp2,nx+2,ny+2);
if (color) disalloc_matrix(tmp3,nx+2,ny+2);

disalloc_matrix(r1, nx+2,ny+2);
disalloc_matrix(r2, nx+2,ny+2);

disalloc_matrix(v1 ,nx+2,ny+2);
/* u2 and u3 for color images only */
if (color) disalloc_matrix(v2 ,nx+2,ny+2);
if (color) disalloc_matrix(v3 ,nx+2,ny+2);
}

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
char  in[80];
char out[80];
char row[80];

int nx,ny;
float hx,hy;

float **u1,**u2,**u3;

int i,j;

float tau,T,t_stop;

unsigned char byte;

int selection;

float sigma,rho;

long pos;

int color;

hx = 1.0;
hy = 1.0;
  
FILE * ptr;

/* get information: color image or grayscale? */
printf("Color [0=No, 1=Yes]:                                    ");
gets(row); sscanf(row,"%d",&color);

/* get input image name */
printf("input image:                                            ");
gets(in);

/* read image data */
if (color==0)
  {
  read_pgm_header(in,&pos,&nx,&ny);
  alloc_matrix(&u1,nx+2,ny+2);
  read_pgm_data(in,pos,u1,nx,ny,1,1);
  }
else
  {
  read_ppm_header(in,&pos,&nx,&ny);
  alloc_matrix(&u1,nx+2,ny+2);
  alloc_matrix(&u2,nx+2,ny+2);
  alloc_matrix(&u3,nx+2,ny+2);
  read_ppm_data_channelwise(in,pos,u1,u2,u3,nx,ny,1,1);
  }

/* get information: what operation should be performed? */  
printf("\n1 = Dilation, 2 = Erosion, 3 = ShockFiltering:        ");
// printf("Select morphological operation:                         ");
gets(row); sscanf(row,"%ld",&selection);

/* check if selection is valid */
if (selection<1||selection>3) {printf("Wrong selection\n"); exit(0);}

/* if shock filter */
if (selection==3)
  {
  /* get sigma */
  printf("Presmoothing scale sigma:                        ");
  gets(row); sscanf(row,"%f",&sigma);

  /* get rho */
  printf("Convolution scale rho:                           ");
  gets(row); sscanf(row,"%f",&rho);
  }

/* get time step size */
printf("Evolution time step:                                    ");
gets(row); sscanf(row,"%f",&tau);

/* get stopping time */
printf("Total evolution time:                                   ");
gets(row); sscanf(row,"%f",&T);

/* get output image name */
printf("output image:                                           ");
gets(out);

t_stop = T;

/* perform selected operation */
switch (selection)
  {
  case 1:
    while (T>tau)
      {
      dilation(u1,nx,ny,hx,hy,tau);
      if (color) dilation(u2,nx,ny,hx,hy,tau);
      if (color) dilation(u3,nx,ny,hx,hy,tau);
      T -= tau;
      }
    dilation(u1,nx,ny,hx,hy,T);
    if (color) dilation(u2,nx,ny,hx,hy,T);
    if (color) dilation(u3,nx,ny,hx,hy,T);
  break;
  case 2:
    while (T>tau)
      {
      erosion(u1,nx,ny,hx,hy,tau);
      if (color) erosion(u2,nx,ny,hx,hy,tau);
      if (color) erosion(u3,nx,ny,hx,hy,tau);
      T -= tau;
      }
    erosion(u1,nx,ny,hx,hy,T);
    if (color) erosion(u2,nx,ny,hx,hy,T);
    if (color) erosion(u3,nx,ny,hx,hy,T);
  break;
  case 3:
    while (T>tau)
      {
      shock_filter(u1,u2,u3,nx,ny,hx,hy,tau,sigma,rho,color);
      T -= tau;
      }
    shock_filter(u1,u2,u3,nx,ny,hx,hy,T,sigma,rho,color);
  break;
  default:
    printf("Wrong selection\n"); exit(0);
  break;
  }    

/* open file and write header (incl. filter parameters) */
ptr = fopen (out, "wb");
fprintf (ptr, "P%d \n",5+color); //P5 for gray, P6 for color images 
fprintf (ptr, "# PDE-based morphology \n");
fprintf (ptr, "# initial image:         %s\n", in);
fprintf (ptr, "# time step:             %f\n",tau);
fprintf (ptr, "# stopping time:         %f\n",t_stop);
fprintf (ptr, "%ld %ld \n255\n", nx, ny);

fclose(ptr);

/* write data */
if (!color)
  write_pgm_data(out,u1,nx,ny,1,1);
else
  write_ppm_data_channelwise(out,u1,u2,u3,nx,ny,1,1);

/* free memory */
disalloc_matrix(u1,nx+2,ny+2);
/* u2 and u3 for color images only */
if (color) disalloc_matrix(u2,nx+2,ny+2);
if (color) disalloc_matrix(u3,nx+2,ny+2);

return(0);
}