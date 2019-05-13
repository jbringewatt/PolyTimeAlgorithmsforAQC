/***********************************
************************************
Author: Jake Bringewatt
Use: Finds the eigenspectrum in the full Hilbert space via exact diagonalization. Takes a standard input file which is for Hamming spherical potentials with constant strength. Can easily be generalized to more complicated potentials.

Dependencies: Lapack, Blas

Compile: gcc -o FullSpaceSolver ./FullSpaceSolver -lm -lblas -llapack
Run: ./FullSpaceSolver wellInputFile paramInputFile

Input: 
Standard well input file (see README) [wellInputFile]:
Has the standard form for all codes:
number of qubits
number of wells
(for each well) 
well location
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


Full space solver parameter input file [paramInputFile]:
Gives relevant parameters for this code:
 
runname - tag for output filenames
include wfs - 0 if no wfs desired, 1 if wfs desired
num energies - number of energies (and wfs if include wfs is 1) to be printed out
sList - 0 if want to telescope, 1 if want to print for a specific list of s, if 1 follow by a list of s values to run problem at
	if 0 next four lines are:  
		numSteps - number of steps in s per range scanned over (recommended to be 10*k+1 for some integer k)
		iMax - maximum number of loops through telescoping procedure
    telescopeParam - 0 for gap within sig=0, 1 for max of 2 possible gaps, 2 for min of 2 possible gaps
    degeneracy - degeneracy of potential function lowest point
	if 1 next lines are list of s values to solve at

e.g
testRun1
1
0
3
21
10

or,

testRun2
0
1
3
0.0
0.1
0.2
0.25
0.30
0.5
0.7
....

Output: 2 csv files:
1. runname.full.energies.output
Gives the energies and gaps at each s where each row is
s, E0 (sigma=0), E1 (sigma=0), E1*(sigma=1), E1-E0, E1*-E0

2. runname.full.wfs.output
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
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int numSteps, int iMax, int telescopeParam,int degeneracy, int numEnergies);
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam, int degeneracy, int numEnergies);
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g,  FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam, int degeneracy, int numEnergies);
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
int binDiff(int a, int b, int B);
int binDiff2(int a, const char* bS, int B);
char *decimalToBinary(int n, int B);



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

  //read in param input file stuff
  getline(&line, &len, paramInputFile);
  char* runname=(char *)malloc((100)*sizeof(char));
  strcpy(runname, line);
  runname[strlen(line)-1] = '\0';
  getline(&line, &len, paramInputFile);
  int getWfs=atoi(line);
  getline(&line, &len, paramInputFile);
  int numEnergies=atoi(line);
  getline(&line, &len, paramInputFile);
  int sList=atoi(line);

  int iMax;
  int numSteps;
  int numS;
  int telescopeParam;
  double* sVals;
  int degeneracy;
  if(sList==0){
	getline(&line, &len, paramInputFile);
  	numSteps=atoi(line);
  	getline(&line, &len, paramInputFile);
  	iMax=atoi(line);
	getline(&line, &len, paramInputFile);
  	telescopeParam=atoi(line);
    getline(&line, &len, paramInputFile);
  	degeneracy=atoi(line);
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
  strcat(energyOutputFilename, ".full.energies.output");  
  energyOutputFile=fopen(energyOutputFilename, "w");
  if(energyOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  }
 
  if(getWfs==1){	  
	char wfOutputFilename[150];
	wfOutputFilename[0]='\0';
	strcat(wfOutputFilename, runname);
	strcat(wfOutputFilename, ".full.wfs.output");
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
   		gsearch(bits, numWells, r, wellLocations, wellWidths, wellDepths, energyOutputFile, wfOutputFile, getWfs, numSteps, iMax, telescopeParam, degeneracy, numEnergies);
	}else{
		//find gap at each s value in sVals
		double g;
		for(i=0; i<numS; i++){
			getStates(sVals[i], bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, degeneracy, numEnergies);
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
Telescopes in on minimum gap, looping until convergence or 10 loops whichever comes first
*/
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int numSteps, int iMax, int telescopeParam, int degeneracy, int numEnergies) {
	double mingap, gapOld;
	int i;
	smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam, degeneracy, numEnergies);	
	printf("range: [%g, %g]\n", r.bottom, r.top);	
    gapOld=1;
    int tol=1e-8;
    int tol2=2e-16;
	int iter=0;
    while(fabs(gapOld-mingap)>tol){ 
		gapOld=mingap;
        if(fabs(r.bottom-r.top)<tol2||iter==iMax) break;
        smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, &mingap, getWfs, numSteps, telescopeParam, degeneracy, numEnergies);
		iter++;       
    }				
}

/*
Performs an even scan of numSteps over a range calculating the gaps, returns next range in r parameter
*/
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths,FILE* energyOutputFile, FILE* wfOutputFile, double* mingap, int getWfs, int numSteps, int telescopeParam, int degeneracy, int numEnergies) {
	int i, mini;
	double s, g;
	double delta;
	range returnval;
	delta = (r->top-r->bottom)/((double)numSteps-1.0);
	getStates(r->bottom, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, degeneracy, numEnergies);
	mini = 0;
	double mingaprange=g;
	for(i = 1; i < numSteps; i++) {
		s = r->bottom + delta*(double)i;
		getStates(s, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, getWfs, telescopeParam, degeneracy, numEnergies);
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
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam, int degeneracy, int numEnergies){
	printf("s: %f\n", s);
    int n = pow(2, bits), lda = n, info, lwork;
    double wkopt;
    double* work;
    /* Local arrays */
    double w[n];
    double* a = malloc(sizeof(double)*lda*n);
    //zero out matrix
    int j;
    for(j=0; j<n*lda; j++){
        a[j]=0.0;
    }
    int i, col, row=-1;
    for(j=0; j<n*lda; j++){
        col=j%n;
        if(col==0) row ++;
        if(col<row){
            if(binDiff(row, col, bits)==1){ 
                a[j]=-1.0*(1-s)/bits;
            }else{
                a[j]=0.0;
            }
        }else if(col==row){	
            for(i=0; i<numWells; i++){
                if(binDiff2(row, wellLocations[i], bits)<wellWidths[i]){
                    a[j]=wellDepths[i]*s;
                    break;
                }else{
                    a[j]=0.0;
                }				
            }				
        }			
    }
            
    /* Executable statements */
    // print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
    /* Query and allocate the optimal workspace */
    lwork = -1;
    dsyev_( "Vectors", "Upper", &n, a, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
    /* Solve eigenproblem */
    dsyev_( "Vectors", "Upper", &n, a, &lda, w, work, &lwork, &info );
    /* Check for convergence */
    if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }
    //write out energies
	fprintf(energyOutputFile, "% 1.15f, ", s);
    for(i=0; i<numEnergies; i++){
        fprintf(energyOutputFile, "% 2.14f, ", w[i]);      
    }
      fprintf(energyOutputFile, "\n");   

    //write out wfs if desired
    if(getWfs==1){
        //write out wfs
		fprintf(wfOutputFile, "% 1.15f\n", s);  
        double fac1, fac2;
		for(j=0; j<numEnergies; j++){
            if(a[n*j]>=-1e-15) fac1=1.0;
			else fac1=-1.0;	
      for(i=0; i<n; i++){                
			    fprintf(wfOutputFile, "% 1.15f, ", fac1*a[j*n+i]);
            }
            fprintf(wfOutputFile, "\n");  		
		}  
    }

    /* Free workspace */
    free(a);
    free( (void*)work );  


    if(telescopeParam==0)
		*g=w[1]-w[0];
	else if(telescopeParam==1)
		*g=max(w[1]-w[0], w[degeneracy]-w[degeneracy-1]);
	else
		*g=min(w[1]-w[0], w[degeneracy]-w[degeneracy-1]);	
}	


/*Returns the difference in binary representation of 2 integers*/
int binDiff(int a, int b, int B){
	char* aS;
	char* bS;
	aS=decimalToBinary(a, B);
	bS=decimalToBinary(b, B);
	int i;
	int numDiff=0;
	for(i=0; i<B; i++){
		if(aS[i]!=bS[i]) numDiff++;
	}
	free(aS);
	free(bS);
	return numDiff;
}

/*Returns the difference in binary representation of an integer and a binary string*/
int binDiff2(int a, const char* bS, int B){
	char *aS;
	aS=decimalToBinary(a, B);
	int i;
	int numDiff=0;
	for(i=0; i<B; i++){
		if(aS[i]!=bS[i]) numDiff++;
	}
	free(aS);
	return numDiff;
}

/* Function to convert a decimal number to binary string,*/
char *decimalToBinary(int n, int B)
{
   int c, d, count;
   char *pointer; 
   count = 0;
   pointer = (char*)malloc(B+1); 
   if ( pointer == NULL )
      exit(EXIT_FAILURE); 
   for ( c = B-1 ; c >= 0 ; c-- )
   {
      d = n >> c; 
      if ( d & 1 )
         *(pointer+count) = 1 + '0';
      else
         *(pointer+count) = 0 + '0'; 
      count++;
   }
   *(pointer+count) = '\0'; 
   return  pointer;
}