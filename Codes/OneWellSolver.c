/***********************************
************************************
Author: Jake Bringewatt
Use: Finds the ground state and first excited state candidate in the sigma = 0 subspace, and one of the n-1 degenerate first excited state candidates in the sigma=1 subspace for a single Hamming symmetric well. Takes a standard input file which is for a Hamming spherical potential with constant strength. Can easily be generalized to more complicated potentials.

Dependencies: Lapack, Blas

Compile: gcc -o OneWellSolver ./OneWellSolver -lm -lblas -llapack
Run: ./OneWellSolver wellInputFile paramInputFile

Input: 
Standard well input file (see README) [wellInputFile]:
Has the standard form for all codes:
number of qubits
number of wells (will return with error if not 1)
well location
well depth (double)(<0)
well width (integer)(in [1,n])

e.g.
10
1
0000011100
-2.1
1

One well solver parameter input file [paramInputFile]:
Gives relevant parameters for this code:
 
runname - tag for output filenames
include wfs - 0 if no wfs desired, 1 if wfs desired
sList - 0 if want to telescope, 1 if want to print for a specific list of s, if 1 follow by a list of s values to run problem at
	if 0 next lines are:  
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
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam);
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g,  FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam);
void getStatesHelper(int width, double depth,int wellNum, int subspace,int bits,int numWells,double s,double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs);
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
double potOneWell(int bits,int width,double depth, int w);
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
  
  if(numWells!=1){
      printf("ERROR: More than one well in input file.\n");
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
  strcat(energyOutputFilename, ".ow.energies.output"); 
	energyOutputFile=fopen(energyOutputFilename, "w");
  if(energyOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  }
 
  if(getWfs==1){	  
	char wfOutputFilename[150];
	wfOutputFilename[0]='\0';
	strcat(wfOutputFilename, runname);
	strcat(wfOutputFilename, ".ow.wfs.output");
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
		double g;
		for(i=0; i<numS; i++){
			getStates(sVals[i], bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam);
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
Shifts basis so well is centered at all zeros bitstring
*/
void basisShift(char** wellLocations, int n){	
	int i;
	for(i=0; i<n; i++){
		if(wellLocations[0][i]=='1'){
			wellLocations[0][i]='0';						
		}		
	}	
}

/*
Telescopes in on minimum gap, looping until convergence or 10 loops whichever comes first
*/
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int numSteps, int iMax, int telescopeParam) {
	double mingap, gapOld;
	int i;
	smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam);	
	printf("range: [%g, %g]\n", r.bottom, r.top);	
    gapOld=1;
    int tol=1e-8;
    int tol2=2e-16;
	int iter=0;
    while(fabs(gapOld-mingap)>tol){ 
		gapOld=mingap;
        if(fabs(r.bottom-r.top)<tol2||iter==iMax) break;
        smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam);
		iter++;       
    }				
}

/*
Performs an even scan of numSteps over a range calculating the gaps, returns next range in r parameter
*/
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths,FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam) {
	int i, mini;
	double s, g;
	double delta;
	range returnval;
	delta = (r->top-r->bottom)/((double)numSteps-1.0);
	getStates(r->bottom, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam);
	mini = 0;
	double mingaprange=g;
	for(i = 1; i < numSteps; i++) {
		s = r->bottom + delta*(double)i;
		getStates(s, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam);
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
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam){
	printf("s: %f\n", s);
    double E0;
    double gaps[2];	
	int j;
    for(j=0; j<2; j++){   //check both subspaces         
	    getStatesHelper(wellWidths[0], wellDepths[0], 0, j, bits, numWells, s, gaps, &E0, energyOutputFile, wfOutputFile, getWfs);
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
void getStatesHelper(int width, double depth,int wellNum, int subspace,int bits,int numWells,double s,double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs){
	int dim=(bits+1-2*subspace);
	int i,j;	
	double* H=(double*)malloc(dim*dim*sizeof(double)); //Full H
	int x, y;
	//fill Operator Hamiltonian
	for(x=0; x<dim; x++){
		for(y=0; y<dim; y++){
			if(x==y){
				H[x*dim+y]=s*potOneWell(bits, width, depth, x+subspace);			
			}
			else if(y==x+1){
				H[x*dim+y]=(1-s)*hopUp(bits, x+subspace, bits, subspace);				
			}		
			else if(y==x-1){
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
    if(subspace==0){
		*E0=w[0];
        fprintf(energyOutputFile, "%1.15f, %1.15f, %1.15f, ", s, w[0], w[1]);
		gaps[0]=w[1]-w[0];       
    }
    else{
		gaps[1]=w[0]-*E0;   
		fprintf(energyOutputFile, "%1.15f, %1.15f, %1.15f\n", w[0], gaps[0], gaps[1]);	
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
        if(subspace==0){
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
Get potential for one well of fixed width
*/
double potOneWell(int bits,int width,double depth, int w){
	if(w<width) return depth;
	return 0.0;
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