#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define PI 3.141592653589793
#define adiiterations 5000

int idum,*iup,*idown,*jup,*jdown,lattice_size_y,lattice_size_x;
float **D,delrho,deltax,deltay,**load,**deflect,**deflect2,te;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
      m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
     m[i] -= ncl;
      }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

void setupgridneighborsperiodic()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=lattice_size_x;
     iup[lattice_size_x]=1;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=lattice_size_y;
     jup[lattice_size_y]=1;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
        unsigned long *v;

        v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
        return v-nl+NR_END;
}

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}

void banbks(a,n,m1,m2,al,indx,b)
float **a,**al,b[];
int m1,m2;
unsigned long indx[],n;
{
        unsigned long i,k,l;
        int mm;
        float dum;

        mm=m1+m2+1;
        l=m1;
        for (k=1;k<=n;k++) {
                i=indx[k];
                if (i != k) SWAP(b[k],b[i])
                if (l < n) l++;
                for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
        }
        l=1;
        for (i=n;i>=1;i--) {
                dum=b[i];
                for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];
                b[i]=dum/a[i][1];
                if (l < mm) l++;
        }
}

#define TINY 1.0e-20

void bandec(a,n,m1,m2,al,indx,d)
float **a,**al,*d;
int m1,m2;
unsigned long indx[],n;
{
        unsigned long i,j,k,l;
        int mm;
        float dum;

        mm=m1+m2+1;
        l=m1;
        for (i=1;i<=m1;i++) {
                for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
                l--;
                for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
        }
        *d=1.0;
        l=m1;
        for (k=1;k<=n;k++) {
                dum=a[k][1];
                i=k;
                if (l < n) l++;
                for (j=k+1;j<=l;j++) {
                        if (fabs(a[j][1]) > fabs(dum)) {
                                dum=a[j][1];
                                i=j;
                        }
                }
                indx[k]=i;
                if (dum == 0.0) a[k][1]=TINY;
                if (i != k) {
                        *d = -(*d);
                        for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
                }
                for (i=k+1;i<=l;i++) {
                        dum=a[i][1]/a[k][1];
                        al[k][i-k]=dum;
                        for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
                        a[i][mm]=0.0;
                }
        }
}
#undef SWAP
#undef TINY


void solveimpl()
{
     int i,j,count;
     float d,**ax,**ay,*xx,*xy,*rx,*ry,**alx,**aly,tot;
     unsigned long *indx,*indy;

     indx=lvector(1,lattice_size_x);
     indy=lvector(1,lattice_size_y);
     ax=matrix(1,lattice_size_x,1,5);
     ay=matrix(1,lattice_size_y,1,5);
     xx=vector(1,lattice_size_x);
     xy=vector(1,lattice_size_y);
     rx=vector(1,lattice_size_x);
     ry=vector(1,lattice_size_y);
     alx=matrix(1,lattice_size_x,1,2);
     aly=matrix(1,lattice_size_y,1,2);
     for (count=1;count<=adiiterations;count++) 
      {for (j=1;j<=lattice_size_y;j++) 
        for (i=1;i<=lattice_size_x;i++)
         deflect2[i][j]=deflect[i][j];       
       tot=0;
       for (j=1;j<=lattice_size_y;j++) {
        for (i=1;i<=lattice_size_x;i++) 
         if ((i>2)&&(i<lattice_size_x-1))
         {rx[i]=-load[i][j];
          ax[i][1]=D[i][j];
          ax[i][2]=-4*D[i][j];
          ax[i][3]=12*D[i][j]+delrho; 
          rx[i]+=-D[i][j]*(deflect2[i][jup[jup[j]]]+deflect2[i][jdown[jdown[j]]]-
           4*(deflect2[i][jup[j]]+deflect2[i][jdown[j]]));
          ax[i][4]=-4*D[i][j];
          ax[i][5]=D[i][j];}
         else 
          {rx[i]=0;
           ax[i][1]=0;
           ax[i][2]=0;
           ax[i][3]=1;
           ax[i][4]=0;
           ax[i][5]=0;} 
        bandec(ax,lattice_size_x,2,2,alx,indx,&d);
        banbks(ax,lattice_size_x,2,2,alx,indx,rx);
        for (i=1;i<=lattice_size_x;i++)
         deflect[i][j]=rx[i];}
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         deflect2[i][j]=deflect[i][j];
       for (i=1;i<=lattice_size_x;i++) {
        for (j=1;j<=lattice_size_y;j++)
         if ((j>2)&&(j<lattice_size_y-1))
          {ry[j]=-load[i][j];
           ay[j][1]=D[i][j];
           ay[j][2]=-4*D[i][j];
           ay[j][3]=12*D[i][j]+delrho;
           ry[j]+=-D[i][j]*(deflect2[iup[iup[i]]][j]+deflect2[idown[idown[i]]][j]-
            4*(deflect2[iup[i]][j]+deflect2[idown[i]][j]));
          ay[j][4]=-4*D[i][j];
          ay[j][5]=D[i][j];}
         else  
          {ry[j]=0;
           ay[j][1]=0;
           ay[j][2]=0;
           ay[j][3]=1;
           ay[j][4]=0;
           ay[j][5]=0;}
        bandec(ay,lattice_size_y,2,2,aly,indy,&d);
        banbks(ay,lattice_size_y,2,2,aly,indy,ry);
        for (j=1;j<=lattice_size_y;j++)
         deflect[i][j]=ry[j];}}
}

main()
{    FILE *fp1,*fp2,*fp3;
     int i,j;

     fp1=fopen("westtelcomb3km","r");
     fp2=fopen("westsnowlinemask3km","r");
     fp3=fopen("westdeflectcombnew","w");
     lattice_size_x=626;
     lattice_size_y=642;
     deltax=3;    /* km */   
     deltay=3;    /* km */
     delrho=6000; /* kg/m^3 */
     deflect=matrix(1,lattice_size_x,1,lattice_size_y);
     deflect2=matrix(1,lattice_size_x,1,lattice_size_y);
     load=matrix(1,lattice_size_x,1,lattice_size_y);
     D=matrix(1,lattice_size_x,1,lattice_size_y);     
     setupgridneighborsperiodic(); 
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp1,"%f",&te);
        if ((i<5)||(i>lattice_size_x-5)||(j<5)||(j>lattice_size_y-5)) te=0;
        D[i][j]=te*te*te*70000000/(deltax*deltax*deltax*deltax)/
         (12*(1-0.25*0.25)); 
        fscanf(fp2,"%f",&load[i][j]);
        load[i][j]*=0.82;  /* rho_c/rho_m */
        if ((i<5)||(i>lattice_size_x-5)||(j<5)||(j>lattice_size_y-5)) load[i][j]=0;
        deflect[i][j]=load[i][j];}
     solveimpl();
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       fprintf(fp3,"%f\n",deflect[i][j]); 
     fclose(fp1);
     fclose(fp2);
     fclose(fp3);
}     
