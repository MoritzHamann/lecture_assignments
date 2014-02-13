#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "libraries.h"

/*----------------------------------------------------------------------------*/

inline float min(float a,float b){return a<b?a:b;}
inline float max(float a,float b){return a>b?a:b;}
inline float sq(float a){return a*a;}

inline float heaviside_function(float z, float epsilon)
{return 0.5*(1 + (2/M_PI)*atan(z/epsilon));}

/*----------------------------------------------------------------------------*/

void idm_updates
(
float **f,
float **v,
float **idm_update,
int nx,
int ny,
int bx,
int by,
float hx,
float hy
)
{
int i,j;

float coeff;
  
float H;
float u_in,u_out;
float A_in,A_out;

float H_eps = 0.1;

float hx_1 = 1.0 / hx;
float hy_1 = 1.0 / hy;
  
u_in  = 0.0; A_in  = 0.0;
u_out = 0.0; A_out = 0.0;

for (i=bx;i<nx+bx;i++)
for (j=by;j<ny+by;j++)
  {
  H = heaviside_function(v[i][j],H_eps);

  /*! TODO !*/
  /*! compute the mean values "u_in" and "u_out" */
  u_in += f[i][j]*H;
  A_in += H;
  u_out += f[i][j]*(1-H);
  A_out += (1-H);
  /*!      !*/
  }

  u_in  /= A_in;
  u_out /= A_out;

for (i=bx;i<nx+bx;i++)
for (j=by;j<ny+by;j++)
  {
  /*! TODO !*/
  /*! compute the coefficient "coeff" of intensity-driven motion */
    coeff = pow(f[i][j]-u_out,2)-pow(f[i][j]-u_in,2);
  /*!      !*/

  if (coeff>0.0) //dilation
    idm_update[i][j] =
                      sqrt(   (   sq( min(v[i][j] - v[i-1][j] , 0 ) )
                                + sq( max(v[i+1][j] - v[i][j] , 0 ) ) ) * hx_1
                            + (   sq( min(v[i][j] - v[i][j-1] , 0 ) )
                                + sq( max(v[i][j+1] - v[i][j] , 0 ) ) ) * hy_1 )
                          * fabs(coeff);
  else //erosion
    idm_update[i][j] =
                    - sqrt(   (   sq( max(v[i][j] - v[i-1][j] , 0 ) )
                                + sq( min(v[i+1][j] - v[i][j] , 0 ) ) ) * hx_1
                            + (   sq( max(v[i][j] - v[i][j-1] , 0 ) )
                                + sq( min(v[i][j+1] - v[i][j] , 0 ) ) ) * hy_1 )
                          * fabs(coeff);
  }
return;
}

/*----------------------------------------------------------------------------*/

void mcm_updates
(
float **v,
float **mcm_update,
int nx,
int ny,
int bx,
int by,
float hx,
float hy
)
{
int i,j;

float vx;
float vy;
float vxx;
float vxy;
float vyy;
float nabla_v_sq;

float hx_1 = 1.0 / hx;
float hx_2 = 1.0 / ( hx * hx );
float hy_1 = 1.0 / hy;
float hy_2 = 1.0 / ( hy * hy );
float hxy4 = 1.0 / ( 4 * hx * hy );

for (i=bx;i<nx+bx;i++)
for (j=by;j<ny+by;j++)
  {
  /*! TODO !*/
  /*! compute the required derivatives vx,vy,vxx,vxy,vyy !*/
    // finite difference from english wikipedia
    vx = (v[i+1][j]-v[i-1][j])*hx_1;
    vy = (v[i][j+1]-v[i][j-1])*hy_1;
    vxx = (v[i+1][j]-2*v[i][j]+v[i-1][j])*hx_2;
    vxy = (v[i+1][j+1]-v[i+1][j-1]-v[i-1][j+1]+v[i-1][j-1])*hxy4;
    vyy = (v[i][j+1]-2*v[i][j]+v[i][j-1])*hy_2;
  /*!      !*/

  nabla_v_sq = vx * vx + vy * vy;

  if (nabla_v_sq<10e-5)
    {
    mcm_update[i][j] = 0.0;
    }
  else
    {
    /*! TODO !*/
    /*! compute the mean curvature motion update term !*/
      mcm_update[i][j] = (pow(vy,2)*vxx-2*vx*vy*vxy + pow(vx,2)*vyy)/(nabla_v_sq);
    /*!      !*/
    }
  }
return;
}

/*----------------------------------------------------------------------------*/

void mean_curvature_motion
(
float **v,
int nx,
int ny,
int bx,
int by,
float hx,
float hy,
float tau
)
{
int i,j;

float **v_new;
float **mcm_update;

alloc_matrix(&mcm_update,nx+2*bx,ny+2*by);

mirror_bounds_2d(v,nx,ny,bx,by);

mcm_updates(v,mcm_update,nx,ny,bx,by,hx,hy);

for (i=bx;i<nx+bx;i++)
for (j=by;j<ny+by;j++)
  {
  v[i][j] = v[i][j] + tau * mcm_update[i][j];
  }

disalloc_matrix(mcm_update,nx+2*bx,ny+2*by);
}

/*----------------------------------------------------------------------------*/

void chan_vese_segmentation
(
float **f,
float **v,
int nx,
int ny,
int bx,
int by,
float hx,
float hy,
float lambda,
float tau
)
{
int i,j;
  
float **v_new;
float **mcm_update;
float **idm_update;

alloc_matrix(&mcm_update,nx+2*bx,ny+2*by);
alloc_matrix(&idm_update,nx+2*bx,ny+2*by);

mirror_bounds_2d(v,nx,ny,bx,by);

idm_updates(f,v,idm_update,nx,ny,bx,by,hx,hy);
mcm_updates(  v,mcm_update,nx,ny,bx,by,hx,hy);

for (i=bx;i<nx+bx;i++)
for (j=by;j<ny+by;j++)
  {
  /*! TODO !*/
  /*! implement the iteration step using "idm_update" and "mcm_update" !*/
  /*! ATTENTION use 1/lambda as idm coefficient for numerical stability !*/
  /*! using lambda as mcm coefficient leads to wrong results !*/

  /*!      !*/
  }
disalloc_matrix(mcm_update,nx+2*bx,ny+2*by);
disalloc_matrix(idm_update,nx+2*bx,ny+2*by);
}

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
char  in[80];
char out[80];
char row[80];
char interm[80];

long pos;

int i,j;

int nx,ny;
int bx,by;
float hx,hy;

float **f;
float **fs;
float **v;
float **v_out;

float T,tau;
float lambda;
float sigma;

int delta;
int c,n;

printf("\n");
printf("Mean Curvature Motion and Chan-Vese Segmentation\n\n");
printf("***********************************************\n\n");
printf("  1/2012 Copyright by Andres Bruhn,            \n");
printf("  Sebastian Volz and Levi Valgaerts            \n\n");
printf("  All rights reserved. Unauthorized usage,     \n");
printf("  copying, hiring, and selling prohibited.     \n\n");
printf("  Send bug reports to                          \n");
printf("  bruhn@vis.uni-stuttgart.de or                \n");
printf("  volz@vis.uni-stuttgart.de                    \n\n");
printf("***********************************************\n\n");

FILE * ptr;

/* get input image name */
printf("\ninput image:                                            ");
gets(in);

printf("\npresmoothing scale sigma:                               ");
gets(row); sscanf(row,"%f",&sigma);

bx = 1;
by = 1;

hx = 1.0;
hy = 1.0;

/* read in image */
read_pgm_header(in,&pos,&nx,&ny);
alloc_matrix(&f,nx+2*bx,ny+2*by);
read_pgm_data(in,pos,f,nx,ny,bx,by);

alloc_matrix(&v,nx+2*bx,ny+2*by);

printf("\ntime step:                                              ");
gets(row); sscanf(row,"%f",&tau);

printf("\nStopping time:                                          ");
gets(row); sscanf(row,"%f",&T);

printf("\nNumber of iterations between writes:                    ");
gets(row); sscanf(row,"%d",&delta);

if (!strcmp(argv[1],"cv"))
  {
  printf("\nMCM weight lambda:                                      ");
  gets(row); sscanf(row,"%f",&lambda);
  }
lambda *= 65536;

printf("\noutput image:                                           ");
gets(out);

if (!strcmp(argv[1],"cv"))
  {
  alloc_matrix(&fs,nx+2*bx,ny+2*by);
  alloc_matrix(&v_out,nx+2*bx,ny+2*by);
    
  init_level_set_circle(v,nx,ny,bx,by,hx,hy);

  presmooth(f,fs,nx,ny,bx,by,sigma);

  c = 0;
  n = 0;

  while (T>tau)
    {
    chan_vese_segmentation(fs,v,nx,ny,bx,by,hx,hy,lambda,tau);
    T -= tau;

    c++;
    if(c==delta)
      {
      level_set_to_image(v, nx, ny, bx, by, v_out);
      n++;
      if (n<10)      sprintf(interm,"0000%d_%s",n,out);
      else if (n<100) sprintf(interm,"000%d_%s",n,out);
      else if (n<1000) sprintf(interm,"00%d_%s",n,out);
      else if (n<10000) sprintf(interm,"0%d_%s",n,out);
      else               sprintf(interm,"%d_%s",n,out);
      ptr = fopen (interm, "wb");
      fprintf (ptr, "P5 \n");
      fprintf (ptr, "%ld %ld \n255\n", nx, ny);
      fclose(ptr);
      write_pgm_data(interm,v_out,nx,ny,bx,by);
      c=0;
      }
    }
  chan_vese_segmentation(fs,v,nx,ny,bx,by,hx,hy,lambda,T);

  level_set_to_image(v, nx, ny, bx, by, v_out);

  sprintf(row,"final-%s",out);
  ptr = fopen (row, "wb");
  fprintf (ptr, "P5 \n");
  fprintf (ptr, "# chan-vese segmentation\n");
  fprintf (ptr, "# initial image:         %s\n", in);
  fprintf (ptr, "# time step:             %f\n",tau);
  fprintf (ptr, "# stopping time:         %f\n",T);
  fprintf (ptr, "# lambda:                %f\n",lambda);
  fprintf (ptr, "%ld %ld \n255\n", nx, ny);
  fclose(ptr);

  write_pgm_data(row,v_out,nx,ny,bx,by);

  mirror_bounds_2d(v_out,nx,ny,bx,by);
  
  for (i=bx;i<nx+bx;i++)
  for (j=bx;j<nx+bx;j++)
    {
    if (v_out[i][j]==255.0)
      {
      if (  (v_out[i-1][j  ]==0)
          ||(v_out[i+1][j  ]==0)
          ||(v_out[i  ][j-1]==0)
          ||(v_out[i  ][j+1]==0)
          ||(v_out[i-1][j-1]==0)
          ||(v_out[i-1][j+1]==0)
          ||(v_out[i+1][j-1]==0)
          ||(v_out[i+1][j+1]==0))
        {
        fs[i][j] = 255.0;
        f[i][j] = 0.0;
        }
      else
        fs[i][j] = f[i][j];
      }        
    }

  sprintf(row,"overlay-%s",out);
  ptr = fopen (row, "wb");
  fprintf (ptr, "P6 \n");
  fprintf (ptr, "# chan-vese segmentation\n");
  fprintf (ptr, "# initial image:         %s\n", in);
  fprintf (ptr, "# time step:             %f\n",tau);
  fprintf (ptr, "# stopping time:         %f\n",T);
  fprintf (ptr, "# lambda:                %f\n",lambda);
  fprintf (ptr, "%ld %ld \n255\n", nx, ny);
  fclose(ptr);
    
  write_ppm_data_channelwise(row,fs,fs,f,nx,ny,bx,by);
  
  disalloc_matrix(fs,nx+2*bx,ny+2*by);
  disalloc_matrix(v_out,nx+2*bx,ny+2*by);
  }
else if (!strcmp(argv[1],"mcm"))
  {
  presmooth(f,v,nx,ny,bx,by,sigma);

  c = 0;
  n = 0;

  while (T>tau)
    {
    mean_curvature_motion(v,nx,ny,bx,by,hx,hy,tau);
    T -= tau;

    c++;
    if(c==delta)
      {
      n++;
      if (n<10)      sprintf(interm,"0000%d_%s",n,out);
      else if (n<100) sprintf(interm,"000%d_%s",n,out);
      else if (n<1000) sprintf(interm,"00%d_%s",n,out);
      else if (n<10000) sprintf(interm,"0%d_%s",n,out);
      else               sprintf(interm,"%d_%s",n,out);
      ptr = fopen (interm, "wb");
      fprintf (ptr, "P5 \n");
      fprintf (ptr, "%ld %ld \n255\n", nx, ny);
      fclose(ptr);
      write_pgm_data(interm,v,nx,ny,bx,by);
      c=0;
      }
    }
  mean_curvature_motion(v,nx,ny,bx,by,hx,hy,T);

  sprintf(row,"final-%s",out);
  ptr = fopen (row, "wb");
  fprintf (ptr, "P5 \n");
  fprintf (ptr, "# mean curvature motion\n");
  fprintf (ptr, "# initial image:         %s\n", in);
  fprintf (ptr, "# time step:             %f\n",tau);
  fprintf (ptr, "# stopping time:         %f\n",T);
  fprintf (ptr, "%ld %ld \n255\n", nx, ny);
  fclose(ptr);

  write_pgm_data(row,v,nx,ny,bx,by);
  }  

disalloc_matrix(f,nx+2*bx,ny+2*by);
disalloc_matrix(v,nx+2*bx,ny+2*by);

return 0;
}
