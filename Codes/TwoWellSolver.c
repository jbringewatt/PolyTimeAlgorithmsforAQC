/***********************************
************************************
Author: Jake Bringewatt
Date:2/11/2019
Use: Finds the ground state and first excited state candidate in the sigma = 0 subspace, and one of the n-1 degenerate first excited state candidates in the sigma=1 (sigma_1=0, sigma_2=1 or sigma_1=1, sigma_2=1) subspace for two Hamming symmetric wells. Takes a standard input file which is for a Hamming spherical potential with constant strength. Can easily be generalized to more complicated potentials.

Dependencies: Lapack, Blas

Compile: gcc -o TwoWellSolver ./TwoWellSolver -lm -lblas -llapack
Run: ./TwoWellSolver wellInputFile paramInputFile

Input: 
Standard well input file (see README) [wellInputFile]:
Has the standard form for all codes:
number of qubits
number of wells (will return with error if not 2)
well depth (double)(<0)
well width (integer)(in [1,n])

e.g.
10
2
1111011101
-2.1
1
0000101110
-2
2

Two well solver parameter input file [paramInputFile]:
Gives relevant parameters for this code:
 
runname - tag for output filenames
include wfs - 0 if no wfs desired, 1 if wfs desired
sList - 0 if want to telescope, 1 if want to print for a specific list of s, if 1 follow by a list of s values to run problem at
	if 0 next two lines are:  
		numSteps - number of steps in s per range scanned over (recommended to be 10*k+1 for some integer k)
		iMax - maximum number of loops through telescoping procedure
        telescopeParam - 0 for gap within sig=0, 1 for max of 2 possible gaps, 2 for min of 2 possible gaps
	if 1 next lines are list of s values to solve at

e.g
testRun1
1
0
21
10

or,

testRun2
0
1
0.0
0.1
0.2
0.25
0.30
0.5
0.7
....

Output: 2 csv files:
1. runname.ow.energies.output
Gives the energies and gaps at each s where each row is
s, E0 (sigma=0), E1 (sigma=0), E1*(sigma=1), E1-E0, E1*-E0

2. runname.ow.wfs.output
For each s gives the following n+1 rows
E0 Wf
E1 Wf
E1* (gamma=0) Wf (same for all gamma)


************************************
***********************************/

//IMPORT STATEMENTS
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<limits.h>
#include<string.h>
#include<float.h>

//Define structures
typedef struct {
  double bottom;
  double top;
}range;

//Define functions
#define max(x,y) ((x) >= (y)) ? (x) : (y)
#define min(x,y) ((x) < (y)) ? (x) : (y)
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int numSteps, int iMax, int telescopeParam);
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam, double* E0);
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g,  FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam, double* E0);
void getStatesHelper(char** wellLocations, int* wellWidths, double* wellDepths,int subspace1, int subspace2, int bits, int numWells,double s, double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs);
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
double pot(int n1, int n2, int h1,int h2,int well1Width,int well2Width, double well1Depth,double well2Depth);
double hopUp(int n,int h,int bits, int sig);
void basisShift(char** wellLocations, int n);


/* Main driver, entry point for program*/
int main(int argc, char *argv[]){
	if(argc != 3){
	  printf("Usage: OneWellSolver wellInputFile paramInputFile\n");
	  return 0;
  }  
  //open input files
  FILE* wellInputFile;
  FILE* paramInputFile;
  wellInputFile=fopen(argv[1], "r");
  paramInputFile=fopen(argv[2], "r");
  if(wellInputFile==NULL || paramInputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  }

  //Read in number of qubits and number of wells
  char *buffer;
  size_t bufsize;
  bufsize = 1;
  buffer = (char *)malloc(bufsize * sizeof(char));
  getline(&buffer, &bufsize, wellInputFile);
  int bits=atoi(buffer);
  getline(&buffer, &bufsize, wellInputFile);
  int numWells=atoi(buffer);
  free(buffer);
  
  if(numWells!=2){
      printf("ERROR: Not two wells in input file.\n");
	  exit(0);
  }

  //get well locations/parameters
  size_t len=1;
  char* line=NULL;
  int i;
  char* wellLocations[numWells];
  double* wellDepths=(double *)malloc(numWells*sizeof(double));
  int* wellWidths=(int *)malloc(numWells*sizeof(int));
  char* wellstring;
  for(i=0; i<numWells; i++){
	getline(&line, &len, wellInputFile);
	wellstring=(char *)malloc((bits+1)*sizeof(char));
	strcpy(wellstring, line);
	wellLocations[i]=wellstring;
	printf("Location: %s", wellLocations[i]);
	getline(&line, &len, wellInputFile);
	wellDepths[i]=atof(line);
	printf("Depth: %f\n", wellDepths[i]);
	getline(&line, &len, wellInputFile);
	wellWidths[i]=atoi(line);
	printf("Width: %i\n", wellWidths[i]);	
  }    
  //shift basis (i.e. put well at all zeros bitstring)
  basisShift(wellLocations, bits);

  //read in param input file stuff
  getline(&line, &len, paramInputFile);
  char* runname=(char *)malloc((100)*sizeof(char));
  strcpy(runname, line);
  runname[strlen(line)-1] = '\0';
  getline(&line, &len, paramInputFile);
  int getWfs=atoi(line);
  getline(&line, &len, paramInputFile);
  int sList=atoi(line);

  int iMax;
  int numSteps;
  int numS;
  int telescopeParam;
  double* sVals;
  if(sList==0){
	getline(&line, &len, paramInputFile);
  	numSteps=atoi(line);
  	getline(&line, &len, paramInputFile);
  	iMax=atoi(line);
	getline(&line, &len, paramInputFile);
  	telescopeParam=atoi(line);
  }else{
	getline(&line, &len, paramInputFile);
  	numS=atoi(line);
	sVals=(double*) malloc(sizeof(double)*numS); 	
	for(i=0; i<numS; i++){
		getline(&line, &len, paramInputFile);
		sVals[i]=atof(line);		
	}
  }   

  //close input files
  fclose(wellInputFile);
  fclose(paramInputFile);

  //create outputfiles
  FILE* energyOutputFile;
  FILE* wfOutputFile;
  char energyOutputFilename[150];
	energyOutputFilename[0]='\0';
  strcat(energyOutputFilename, runname);
  strcat(energyOutputFilename, ".tw.energies.output");  
  energyOutputFile=fopen(energyOutputFilename, "w");
  if(energyOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  }
 
  if(getWfs==1){	  
	char wfOutputFilename[150];
	wfOutputFilename[0]='\0';
	strcat(wfOutputFilename, runname);
	strcat(wfOutputFilename, ".tw.wfs.output");
	wfOutputFile=fopen(wfOutputFilename, "w");
	if(wfOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  	}
  }

	if(sList==0){
		//find gap as a function of s by telescoping in on region of minimum gap
		range r;
  		r.top=1.0;
   		r.bottom=0.0;
   		gsearch(bits, numWells, r, wellLocations, wellWidths, wellDepths, energyOutputFile, wfOutputFile, getWfs, numSteps, iMax, telescopeParam);
	}else{
		//find gap at each s value in sVals
		double g, E0;
		for(i=0; i<numS; i++){
			getStates(sVals[i], bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, &E0);
		}
	}   
      
  //close output files
  if(getWfs==1)
  	fclose(wfOutputFile);
  fclose(energyOutputFile);

  //free resources
  for(i=0; i<numWells; i++){
	  free(wellLocations[i]);
  }
  free(wellDepths);
  free(wellWidths);
  free(line);
  free(runname);
  if(sList==1)
  	free(sVals);
  exit( 0 );
}

/*
Shifts basis so one well is on the all zeros bitstring 
*/
void basisShift(char** wellLocations, int n){	
	int i;
	for(i=0; i<n; i++){
		if(wellLocations[0][i]=='1'&&wellLocations[1][i]=='1'){
			wellLocations[0][i]='0';
			wellLocations[1][i]='0';			
		}
		else if(wellLocations[0][i]=='1'&&wellLocations[1][i]=='0'){
			wellLocations[0][i]='0';
			wellLocations[1][i]='1';			
		}		
	}	
}

/*
Telescopes in on minimum gap, looping until convergence or 10 loops whichever comes first
*/
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int numSteps, int iMax, int telescopeParam) {
	double mingap, gapOld, E0, E0Old;
	int i;
	smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam, &E0);	
	printf("range: [%g, %g]\n", r.bottom, r.top);	
    gapOld=1;
		E0Old=10;
    int tol=2e-16;
    int tol2=2e-16;
	int iter=0;
    while(fabs(E0Old-E0)>tol){ 
			gapOld=mingap;
			E0Old=E0;
      if(fabs(r.bottom-r.top)<tol2||iter==iMax) break;
      smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam, &E0);
			iter++;       
    }				
}

/*
Performs an even scan of numSteps over a range calculating the gaps, returns next range in r parameter
*/
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths,FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam, double* E0) {
	int i, mini;
	double s, g;
	double delta;
	range returnval;
	delta = (r->top-r->bottom)/((double)numSteps-1.0);
	getStates(r->bottom, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, E0);
	mini = 0;
	double mingaprange=g;
	for(i = 1; i < numSteps; i++) {
		s = r->bottom + delta*(double)i;
		getStates(s, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, E0);
		if(g < mingaprange) {
		  mingaprange = g;
		  mini = i;
		}
	}
	*mingap=mingaprange;
	if(mini == numSteps-1) {
		r->top = r->top;
		r->bottom = r->top - delta;
		return;
	}
	if(mini == 0) {
		r->bottom = r->bottom;
		r->top = r->bottom + delta;	
		return;	
	}
    r->top = r->bottom+delta*(double)(mini+1);
	r->bottom = r->bottom+delta*(double)(mini-1);		
}

/*
Calculates the energies and if desired wfs as a given s for sigma = 0 and sigma =1 subspaces
*/
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam, double* E0){
	printf("s: %f\n", s);
  
    double gaps[4];	
	int i,j;
    for(i=0; i<2; i++){
        for(j=0; j<2; j++){   //check both subspaces for both sig1 and sig2   
	        getStatesHelper(wellLocations, wellWidths, wellDepths, i, j, bits, numWells, s, gaps, E0, energyOutputFile, wfOutputFile, getWfs);
        }
    }
    	
	if(telescopeParam==0)
		*g=gaps[0];
	else if(telescopeParam==1)
		*g=max(gaps[0], gaps[1]);
	else
		*g=min(gaps[0], gaps[1]);		
}	


/*
Used to find energies and wfs for a given subspace
*/
void getStatesHelper(char** wellLocations, int* wellWidths, double* wellDepths,int subspace1, int subspace2, int bits, int numWells,double s, double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs){
	//find well separations
	int i, j, k, n1, n2;
	//Find n1 and get list of n1 indicies
	n1=0;
	n2=0;
	for(k=0; k<bits; k++){		
		if(wellLocations[1][k]=='1'){
			n1++;
		}else{
			n2++;
		}		
	}
	int sig1=subspace1, sig2=subspace2;
	if(n1==0||n1==1) sig1=0;
	if(n2==0||n2==1) sig2=0;
	int dim=(n1+1-2*sig1)*(n2+1-2*sig2);
	double well1Depth=wellDepths[0];
	double well2Depth=wellDepths[1];
	int well1Width=wellWidths[0];
	int well2Width=wellWidths[1];
	double* H=(double *)malloc(dim*dim*sizeof(double));
	int h1=0, h2=0;
	int x, y;
	//fill Operator Hamiltonian
	for(x=0; x<dim; x++){
		for(y=0; y<dim; y++){
			h1=x/(n2+1-2*sig2)+sig1;
			h2=x%(n2+1-2*sig2)+sig2;	
			//printf("%i, %i, %i, %i\n", sig1, sig2, h1, h2);	
			if(x==y){
				H[x*dim+y]=s*pot(n1, n2, h1, h2, well1Width, well2Width, well1Depth, well2Depth);					
			}
			else if(y==x+1&&n2!=0){
				H[x*dim+y]=(1-s)*hopUp(n2, h2, bits, sig2);				
			}		
			else if(y==x-1&&n2!=0){
				H[x*dim+y]=H[y*dim+x];			
			}
			else if(y==x+n2+1-2*sig2&&n1!=0){
				H[x*dim+y]=(1-s)*hopUp(n1, h1, bits, sig1);				
			}
			else if(y==x-n2-1+2*sig2&&n1!=0){
				H[x*dim+y]=H[y*dim+x];
			}
			else{
				H[x*dim+y]=0.0;				
			}			
		}		
	}
	int n = dim, lda = n, info, lwork;
	double wkopt;
	double* work1;
	/* Local arrays */
	double w[n];
	/* Query and allocate the optimal workspace */
	lwork = -1;
	char mode;
	if(getWfs==1)
		mode='V';
	else
		mode='N';
	dsyev_(&mode, "Upper", &n, H, &lda, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work1 = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_(&mode, "Upper", &n, H, &lda, w, work1, &lwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
	}	
    if(subspace1==0&&subspace2==0){
		*E0=w[0];
        fprintf(energyOutputFile, "%1.15f, %1.15f, %1.15f, ", s, w[0], w[1]);
		gaps[0]=w[1]-w[0];       
    }
    else if(subspace1==0&&subspace2==1){
        gaps[1]=w[0]-*E0;   
		fprintf(energyOutputFile, "%1.15f, ", w[0]);	
    }
    else if(subspace1==1&&subspace2==0){
        gaps[2]=w[0]-*E0;   
		fprintf(energyOutputFile, "%1.15f, ", w[0]);	
    }
    else{
        gaps[3]=w[0]-*E0;   
		fprintf(energyOutputFile, "%1.15f, %1.15f, %1.15f, %1.15f, %1.15f\n", w[0], gaps[0], gaps[1], gaps[2], gaps[3]);	
    }
    
    if(getWfs==1){
        //write out wfs       
		fprintf(wfOutputFile, "%1.15f\n", s);  
        double fac1, fac2;
		if(H[0]>=-1e-15) fac1=1.0;
		else fac1=-1.0;		
		for(j=0; j<dim; j++){
			fprintf(wfOutputFile, "% 1.15f, ", fac1*H[j]);		
		}            
		fprintf(wfOutputFile, "\n");         
        if(subspace1==0&&subspace2==0){
            double fac1, fac2;
            if(H[dim]>=-1e-15) fac1=1.0;
            else fac1=-1.0;		
            for(j=0; j<n; j++){
                fprintf(wfOutputFile, "% 1.15f, ", fac1*H[dim+j]);
            }
			fprintf(wfOutputFile, "\n"); 
        }         
    }       
	free(work1);
	free(H);	
}

/*
2 well potential
*/
double pot(int n1, int n2, int h1,int h2,int well1Width,int well2Width, double well1Depth,double well2Depth){
	//note must guarantee wells don't overlap beforehand
	//printf("Well 3 Depth: %g\n", well3Depth);
	if(h1+h2<well1Width){
		return well1Depth;
	}
	else if(((n1-h1)+h2)<well2Width){
		return well2Depth;
	}
	return 0;	
}

/*
Hop up in subspace given by sig 
*/
double hopUp(int n,int h,int bits, int sig){
	double dn, dh, dbits, ds;
	dn = (double) n;
	dh = (double) h;
	dbits = (double) bits;
	ds = (double) sig;
	return -1*sqrt((dh-ds+1.0)*(dn-ds-dh))/dbits;
}