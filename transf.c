
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#define PI          3.1415926535897932 
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563)



double *mat(int n, int m)
{
    double *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
        printf ("Error al alocar memoria");
        }
    return p;
}
void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);
    
    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}

double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
void xyz2enu(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}
void ecef2enu(const double *pos, const double *r, double *e)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,r,0.0,e);
}

void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}

static int write_results(const double *xstate, int n, FILE *file){
	int cont;
	for (cont=0;cont <n;cont++){
		if(cont == (n-1)){
			fprintf (file,"%lf\n",xstate[cont]);
		}
		else{	
			
			fprintf (file, "%lf,",xstate[cont]);	
		}
		
	
	}

	return 0;
}

int main(){
	FILE *valores;
	FILE *referenceenu;
	valores=fopen("xyzre.txt","r");
	int counti, i,j;
	double *xp;
	xp= mat(3,1);
	double rbase[3], rrprint[3], posprint[3],enuprint[3];
	rbase[0]=-1641945.9160;
	rbase[1]=-3664809.3072;
	rbase[2]=4940010.6554;
	for (i=0;i<20497;i++){
		for(j=0;j<3;j++){
			fscanf(valores,"%lf",&xp[j]);
		
		}
		for (counti=0;counti<3;counti++) rrprint[counti]= xp[counti]-rbase[counti];
		referenceenu = fopen("referenciasenu.csv","a");
		ecef2pos(rbase,posprint);
		ecef2enu(posprint,rrprint,enuprint);
		write_results(enuprint,3,referenceenu);
		fclose(referenceenu);
	
	}

	fclose(valores);

	return 0;
}





