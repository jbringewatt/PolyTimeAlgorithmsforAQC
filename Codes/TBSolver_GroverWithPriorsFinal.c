/***********************************
************************************
Author: Jake Bringewatt
Use: A tight binding code for individual Hamming spherical wells with fixed depth and width, as given by the well input file.  Wells except the first well listed in this file will be given a schedule of b(s)=1-s instead of b(s)=s. Used for doing a Grover with priors problem as in example 1 of the paper.

Dependencies: Lapack, Blas

Compile: 
gcc -c TBSolver.c
gfortran -c dsygvic.f90 dlaasrt.f90 -llapack -lblas
gfortran -o TBSolver TBSolver.o dsygvic.o dlaasrt.f90 -lm -lblas -llapack
rm *.o
(Run TBCompile.bash)
Run: ./TBSolver wellInputFile paramInputFile

Input: 
Standard well input file (see README) [wellInputFile]:
Has the standard form for all codes:
number of qubits
number of wells 
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

One well solver parameter input file [paramInputFile]:
Gives relevant parameters for this code:
 
runname - tag for output filenames
include wfs - 0 if no wfs desired, 1 if wfs desired
include tbwfs - 0 if individual well wfs not desired, 1 if they are
tbOrder - 0 if zeroth order, 1 if 1st order
etol - tolerance for discarding eigenvalues during generalized eigensolver, if -1 automatically find etol (recommended starting value 1e-3)
sList - 0 if want to telescope, 1 if want to print for a specific list of s, if 1 follow by a list of s values to run problem at
	if 0 next four lines are:  
		numSteps - number of steps in s per range scanned over (recommended to be 10*k+1 for some integer k)
		iMax - maximum number of loops through telescoping procedure
        telescopeParam - 0 for gap within sig=0, 1 for max of 2 possible gaps, 2 for min of 2 possible gaps
		degeneracy - degeneracy of potential function lowest point
	if 1 next lines are number of s values and then list of s values to solve at


e.g
testRun1
1
1
0
1e-3
0
21
10
0
1

or,

testRun2
0
0
1
1e-3
1
0.0
0.1
0.2
0.25
0.30
0.5
0.7
....

Output: 3 csv files:
1. runname.tb.energies.output
Gives the energies and gaps at each s where each row is
s, E0, E1, ... E_n, average Energy Difference, average pertubation ratio, total energy difference, total pertubation ratio (THIS IS THE GOOD ERROR ESTIMATE), average of off diagonal elements, sum of off diagonal elements

2. runname.tb.wfs.output
For each s gives the following 2 rows (from H^(TB))
E0 Wf
E1 Wf

3. runname.tb.tbwfs.output
s, bits, n1, n2,dim, dim1, dim2, wellNum1, wfNum1, wellNum2, wfNum2
wf1
wf2


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
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfsOutputFile, int getWfs, int getTbWfs, int numSteps, int iMax, int telescopeParam, int tbOrder, double* etol, int degeneracy);
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfOutputFile,  double* mingap, int getWfs, int getTbWfs, int numSteps, int telescopeParam, int tbOrder, double* etol, int degeneracy, double* E0);
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g,  FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfOutputFile, int getWfs, int getTbWfs, int telescopeParam, int tbOrder, double* etol, int degeneracy, double* gsE0);
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );
extern void dsygvic_(int* itype, char* jobz, char* uplo, int* n, double* a, int* lda, double* b, int* ldb, double* etol, int* k, double* w, double* work, int* ldwork, double* work2, int* lwork2, int* iwork, int* info);
double pot(int n1, int n2, int h1,int h2,int well1Width,int well2Width, double well1Depth,double well2Depth);
double hopUp(int n,int h,int bits, int sig);
void getInternalGap(int wellWidth, double wellDepth, int bits, int wellNum, int subspace, double s, double* internalGaps, int* wellSubspaces, double* E0, double* E1);
double potOneWell(int bits,int width,double depth, int w);
void generateArrayElements(double s,int bits,int wellNum1,int wellNum2,int* R,int numWells, int* wellWidths,double* wellDepths,int* wellSubspaces,double** psiPsiArray,double** psiHPsiArray, FILE* tbWfOutputFile,int getTbWfs,int tbOrder, int numWfs, int* gammaArray, int* divisionIdx, double** errorMatrix);	
//void generateArrayElements(double s,int bits,int wellNum1,int wellNum2,int wfNum1,int wfNum2,int* R,int numWells,int* wellWidths,double* wellDepths, int* wellSubspaces, double* psiPsi,double* psiHPsi, FILE* tbWfOutputFile, int getTbWfs);
double potAvg(int h1,int h2, int n1, int bits, int* wellWidths,double* wellDepths, int i, int j, int numWells, int* R, double s);
double potAvgErr(int h1,int h2, int n1, int bits, int* wellWidths,double* wellDepths, int i, int j, int numWells, int* R, double s);
int N(int R12, int R23, int n1, int h1, int h2, int w, int bits);
double binomial(int n, int k);
double Gamma(int bits,int sig);
double getError(double* errorMatrix, int dim);

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
  
  if(numWells<=0){
      printf("ERROR: Less than one well in input file.\n");
	  exit(0);
  }

  //get well locations/parameters
  size_t len=bits+1;
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
  int getTbWfs=atoi(line);
  getline(&line, &len, paramInputFile);
  int tbOrder=atoi(line);
  getline(&line, &len, paramInputFile);
  double etol=atof(line);
  getline(&line, &len, paramInputFile);
  int sList=atoi(line);
  

  int iMax;
  int numSteps;
  int numS;
  int telescopeParam;
  int degeneracy;
  double* sVals;
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
  FILE* tbWfOutputFile;
  char energyOutputFilename[150];
  energyOutputFilename[0]='\0';
  strcat(energyOutputFilename, runname);
  strcat(energyOutputFilename, ".tb.energies.output");  
  energyOutputFile=fopen(energyOutputFilename, "w");
  if(energyOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  }
 
  if(getWfs==1){	  
	char wfOutputFilename[150];
	wfOutputFilename[0]='\0';
	strcat(wfOutputFilename, runname);
	strcat(wfOutputFilename, ".tb.wfs.output");
	wfOutputFile=fopen(wfOutputFilename, "w");
	if(wfOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  	}
  }
  if(getTbWfs==1){	  
	char tbWfOutputFilename[150];
	tbWfOutputFilename[0]='\0';
	strcat(tbWfOutputFilename, runname);
	strcat(tbWfOutputFilename, ".tb.tbwfs.output");
	tbWfOutputFile=fopen(tbWfOutputFilename, "w");
	if(tbWfOutputFile==NULL){
	  printf("ERROR: Cannot open file \n");
	  exit(0);
  	}
  }

	if(sList==0){
		//find gap as a function of s by telescoping in on region of minimum gap
		range r;
  		r.top=1.0;
   		r.bottom=0.0;
   		gsearch(bits, numWells, r, wellLocations, wellWidths, wellDepths, energyOutputFile, wfOutputFile, tbWfOutputFile, getWfs, getTbWfs, numSteps, iMax, telescopeParam, tbOrder, &etol, degeneracy);
	}else{
		//find gap at each s value in sVals
		double g, E0;
		for(i=0; i<numS; i++){
			getStates(sVals[i], bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, tbWfOutputFile, getWfs, getTbWfs, telescopeParam, tbOrder, &etol, degeneracy, &E0);
		}
	}   
      
  //close output files
  if(getWfs==1)
  	fclose(wfOutputFile);
  if(getTbWfs==1)
	fclose(tbWfOutputFile);
  fclose(energyOutputFile);

  //free resources
  free(wellstring);
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
void gsearch(int bits, int numWells, range r,char** wellLocations,int* wellWidths,double* wellDepths, FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfOutputFile, int getWfs, int getTbWfs, int numSteps, int iMax, int telescopeParam, int tbOrder, double* etol, int degeneracy) {
	double mingap, gapOld, E0, E0Old;
	int i;
	smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, tbWfOutputFile, &mingap, getWfs, getTbWfs, numSteps, telescopeParam, tbOrder, etol, degeneracy, &E0);	
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
        smin(bits, numWells, &r, wellLocations,wellWidths,wellDepths, energyOutputFile, wfOutputFile, tbWfOutputFile, &mingap, getWfs, getTbWfs, numSteps, telescopeParam, tbOrder, etol, degeneracy, &E0);
		iter++;       
    }				
}


/*
Performs an even scan of numSteps over a range calculating the gaps, returns next range in r parameter
*/
void smin(int bits, int numWells, range* r, char** wellLocations,int* wellWidths,double* wellDepths,FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfOutputFile, double* mingap, int getWfs, int getTbWfs, int numSteps, int telescopeParam, int tbOrder, double* etol, int degeneracy, double* E0) {
	int i, mini;
	double s, g;
	double delta;
	range returnval;
	printf("nSteps: %i\n", numSteps);
	delta = (r->top-r->bottom)/((double)numSteps-1.0);
	getStates(r->bottom, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, tbWfOutputFile, getWfs, getTbWfs, telescopeParam, tbOrder, etol, degeneracy, E0);
	mini = 0;
	double mingaprange=g;
	//printf("mingaprange: %g\n", mingaprange);
	for(i = 1; i < numSteps; i++) {
		s = r->bottom + delta*(double)i;
		getStates(s, bits, numWells, wellLocations,wellWidths,wellDepths, &g, energyOutputFile, wfOutputFile, tbWfOutputFile, getWfs, getTbWfs, telescopeParam, tbOrder, etol, degeneracy, E0);
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
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g, FILE* energyOutputFile, FILE* wfOutputFile, FILE* tbWfOutputFile, int getWfs, int getTbWfs, int telescopeParam, int tbOrder, double* etol, int degeneracy, double* tbE0){
	printf("s: %f\n", s);
    //find well separations
	int *R=(int *)malloc(numWells*numWells*sizeof(int));
	int i, j, k, numDiff;
	int minSep=bits;
	for(i=0; i<numWells; i++){ 
		for(j=0; j<numWells; j++){
			if(i==j){
				R[i*numWells+j]=0;
			}
			else if(j<i){
				R[i*numWells+j]=R[j*numWells+i];
			}
			else{
				numDiff=0;
				for(k=0; k<bits; k++){
                   	if(wellLocations[i][k]!=wellLocations[j][k])
                        numDiff++;	
                }
				R[i*numWells+j]=numDiff;
				if(numDiff<minSep) minSep=numDiff;
			}
		}
	} 	

	//view R
	// for(i=0; i<numWells; i++){
	// 	for(j=0; j<numWells; j++){
	// 		printf("%i ", R[i*numWells+j]);
	// 	}
	// 	printf("\n");
	// }
    
    //get number of wfs
    int numWfsXD=0; //numwfs excluding degeneracies
    if(tbOrder==0){
        numWfsXD=numWells;		
	}
    else if(tbOrder==1){		
		//depends on how many excited states there are
        numWfsXD=numWells*2;
    }
    double E0, E1;
    double internalGaps[numWells];	
	double basisStateEnergies[numWfsXD];
    int wellSubspaces[numWells];
    //Check what internal gaps are even if we don't use first excited states in tight binding, keep track of which excited states/subspaces we need
	for(i=0; i<numWells; i++){
        for(j=0; j<2; j++){
            getInternalGap(wellWidths[i], wellDepths[i], bits, i, j, s, internalGaps, wellSubspaces, &E0, &E1);
        }
		if(tbOrder==0){
			basisStateEnergies[i]=E0;			
		}else{
			basisStateEnergies[i]=E0;			
			basisStateEnergies[numWells+i]=E1;
		}
    }
			
    
	//based on subspace calculation calculate numWfs
	int numWfs=0; //numwfs with degeneracies
	int gammaArray[numWells];
	int divisionIdx[numWells+1];
	divisionIdx[0]=0;
	if(tbOrder==0){
		numWfs=numWells;
		for(i=0; i<numWells; i++){
			divisionIdx[i+1]=divisionIdx[i]+1;
			gammaArray[i]=1;
		}
	}else{
		for(i=0; i<numWells; i++){			
			gammaArray[i]=(int) Gamma(bits, wellSubspaces[i]);
			divisionIdx[i+1]=divisionIdx[i]+gammaArray[i]+1;
			numWfs+=gammaArray[i]+1;
		}
	}
	
	
	double basisStateEnergiesFull[numWfs];
	for(i=0; i<numWells; i++){
		basisStateEnergiesFull[divisionIdx[i]]=basisStateEnergies[i];
		for(j=1; j<=gammaArray[i]; j++){
			basisStateEnergiesFull[divisionIdx[i]+j]=basisStateEnergies[numWells+i];
		}		
	}

	double internalGapsFull[numWfs];
	for(i=0; i<numWells; i++){
		for(j=0; j<=gammaArray[i]; j++){
			internalGapsFull[divisionIdx[i]+j]=internalGaps[i];
		}		
	}
	

	//create matrices H^(TB) and S
	double* psiPsiArray=(double*) malloc(sizeof(double)*numWfs*numWfs);
	double* psiHPsiArray=(double*) malloc(sizeof(double)*numWfs*numWfs);
	double* errorMatrix=(double*) malloc(sizeof(double)*numWfs*numWfs);
	//initialize arrays to all zeros
	for(i=0; i<numWfs; i++){
		for(j=0; j<numWfs; j++){
			psiPsiArray[i*numWfs+j]=0.0;
			psiHPsiArray[i*numWfs+j]=0.0;
		}
	}
	for(i=0; i<numWells; i++){
		for(j=i; j<numWells; j++){
			generateArrayElements(s, bits, i, j, R, numWells, wellWidths, wellDepths, wellSubspaces, &psiPsiArray, &psiHPsiArray, tbWfOutputFile, getTbWfs, tbOrder, numWfs, gammaArray, divisionIdx, &errorMatrix);			
		}				
	}

    //now solve generalized eigenproblem
    //eliminate elements that have strong overlaps
	int flaggedWfs[numWfs];
	int numFlagged=0;
	for(i=0; i<numWfs; i++){
		for(j=0; j<numWfs; j++){
			if(j >= i){
				if(psiPsiArray[i*numWfs+j]>1 & j!=i){
					flaggedWfs[i]=1;
					numFlagged++;
					break;
				}else{
					flaggedWfs[i]=0;
				}
			}
		}
	}
	for(i=0; i<numWells; i++){
		if(wellWidths[i]==1 && tbOrder==1){
			flaggedWfs[divisionIdx[i]+1]=1;
			numFlagged++;
		}
	}

	int numWfsNew=numWfs-numFlagged;
	double* psiPsiArrayN=(double*) malloc(sizeof(double)*numWfsNew*numWfsNew);
	double* psiHPsiArrayN=(double*) malloc(sizeof(double)*numWfsNew*numWfsNew);
	double* errorMatrixN=(double*) malloc(sizeof(double)*numWfsNew*numWfsNew);
	
	int iNew=0;
	int jNew=0;
	for(i=0; i<numWfs; i++){
		if(flaggedWfs[i]==0){
			jNew=0;
			for(j=0; j<numWfs; j++){
				if(flaggedWfs[j]==0){
					psiPsiArrayN[iNew*numWfsNew+jNew]=psiPsiArray[i*numWfs+j];
					psiHPsiArrayN[iNew*numWfsNew+jNew]=psiHPsiArray[i*numWfs+j];
					errorMatrixN[iNew*numWfsNew+jNew]=errorMatrix[i*numWfs+j];
					jNew++;
				}
			}
			iNew++;
		}		
	}

	
	//view TB Hamiltonian
	printf("HTB\n");
	for(i=0; i<numWfsNew; i++){
		for(j=0; j<numWfsNew;j++){
			printf("% 3.5g ", psiHPsiArrayN[i*numWfsNew+j]);
		}
		printf("\n");
	}

	printf("S\n");
	for(i=0; i<numWfsNew; i++){
		for(j=0; j<numWfsNew;j++){
			printf("% 3.5g ", psiPsiArrayN[i*numWfsNew+j]);
		}
		printf("\n");
	}
	
	// //get first order energy pertubation for each wf for error estimates
	double pertDiff[numWfs];
	double pertRatio[numWfs];
	double sum=0, sum1=0;

	for(i=0; i<numWfs; i++){
		if(flaggedWfs[i]==0){
			pertDiff[i]=psiHPsiArray[i*numWfs+i]-basisStateEnergiesFull[i];
			sum+=pertDiff[i]*pertDiff[i];
			pertRatio[i]=pertDiff[i]*pertDiff[i]/internalGapsFull[i];
			sum1+=pertRatio[i]*pertRatio[i];
		}else{
			pertDiff[i]=NAN;
			pertRatio[i]=NAN;
		}		
	}
	sum1=sqrt(sum1);
	double totEnergyDiff=sqrt(sum);
	double totPertRatio=sum1;
	double avgEnergyDiff=sum/numWfsNew;
	double avgPertRatio=sum1/numWfsNew;

	//get off diagonal error estimate
	sum=0.0;
	for(i=0; i<numWfsNew; i++){
		for(j=0; j<numWfsNew; j++){
			if(i!=j){
				sum+=psiPsiArrayN[i*numWfsNew+j]*psiPsiArrayN[i*numWfsNew+j];
			}
		}
	}
	double offDiagErrorEst=sqrt(sum);
	double avgOffDiagErrorEst=sum/numWfsNew;


    int numWfsOld=numWfs;
    numWfs=numWfsNew;
	double workoutput;
    int n=numWfs;
    int lda=numWfs;
	int ldb=numWfs;
   	int itype=1; 
	int lwork=-1;
    int info;
    double w1[n];
    double* work2;
    int iwork[n];
    int kval;
    int ldwork=numWfs;
    double* work=(double*)malloc(sizeof(double)*numWfs*numWfs);

	char mode;
	dsygvic_(&itype, "V", "U", &n, psiHPsiArrayN, &lda, psiPsiArrayN, &ldb, etol, &kval, w1, work, &n, &workoutput, &lwork, iwork, &info);
    
    lwork = (int)workoutput;
    work2 = (double*)malloc( lwork*sizeof(double) );
    dsygvic_(&itype, "V", "U", &n, psiHPsiArrayN, &lda, psiPsiArrayN, &ldb, etol, &kval, w1, work, &n, work2, &lwork, iwork, &info);

	printf("kval: %i\n", kval);

	//write out energies
	fprintf(energyOutputFile, "% 1.15f, ", s);
    for(i=0; i<numWfsOld; i++){
        if(i<kval)
            fprintf(energyOutputFile, "% 2.14f, ", w1[i]);
        else
            fprintf(energyOutputFile, "              nan, ");
    } 

	double h, h1, h2;
	double hX1[numWfs];
	double hX2[numWfs];
	for(i=0; i<numWfs; i++){
		h1=0.0;
		h2=0.0;
		for(j=0; j<numWfs; j++){
			h1+=errorMatrix[i*numWfs+j]*psiHPsiArrayN[j];
			h2+=errorMatrix[i*numWfs+j]*psiHPsiArrayN[numWfs+j];
		}
		hX1[i]=h1;
		hX2[i]=h2;
	}
	h1=0;
	h2=0;
	for(i=0; i<numWfs; i++){
		h1+=psiHPsiArrayN[i]*hX1[i];
		h2+=psiHPsiArrayN[numWfs+i]*hX2[i];
	}
	h=max(h1,h2);
	double evalDiff=w1[2]-w1[1];

	//diagonalize error matrix to get error estimate
	// double errorEst=fabs(getError(errorMatrix, numWfs));
	// if(evalDiff>2*errorEst-h){
	// 	errorEst=h+errorEst*errorEst/(evalDiff+h-errorEst);
	// 	fprintf(energyOutputFile, "% 2.14f, ", errorEst);  
	// }else{		
	// 	fprintf(energyOutputFile, "% 2.14f, ", errorEst);  
	// }
	// //fprintf(energyOutputFile, "% 2.14f\n", errorEst);

	// fprintf(energyOutputFile, "% 2.14f, ", offDiagErrorEst);


	// //print out pertubation theory error estimates per well
	// for(i=0; i<numWfsOld; i++){
	// 	if(isnan(pertDiff[i])){
	// 		fprintf(energyOutputFile, "              nan,               nan, ");
	// 	}else if(isnan(pertRatio[i])){
	// 		fprintf(energyOutputFile, "% 2.14f,               nan, ", pertDiff[i]);
	// 	}
	// 	else{				
	// 		fprintf(energyOutputFile, "% 2.14f, % 2.14f, ", pertDiff[i], pertRatio[i]);
	// 	}
	// }

	// //print out internal gap energies
	// sum=0.0;
	//  for(i=0; i<numWfsOld; i++){
    //     fprintf(energyOutputFile, "% 2.14f, ", internalGaps[i]);  
	// 	sum+=internalGaps[i];      
    // }
	// //and the average
	// fprintf(energyOutputFile, "% 2.14f, ", sum/numWfsOld);  

	fprintf(energyOutputFile, "% 2.14f, % 2.14f, % 2.14f, % 2.14f, % 2.14f, % 2.14f \n", avgEnergyDiff, avgPertRatio, totEnergyDiff, totPertRatio, avgOffDiagErrorEst, offDiagErrorEst);
		
	if(getWfs==1){
		//write out wfs
		fprintf(wfOutputFile, "% 1.15f\n", s);  
        double fac1, fac2;			
		for(j=0; j<kval; j++){
			if(psiHPsiArrayN[kval*j]>=-1e-15) fac1=1.0;
			else fac1=-1.0;	
			for(i=0; i<kval; i++){				
				fprintf(wfOutputFile, "% 1.15f, ", fac1*psiHPsiArrayN[j*kval+i]);
			}				
			fprintf(wfOutputFile, "\n"); 	
		}  
	}
	
	*tbE0=w1[0];
	if(telescopeParam==0)
		*g=w1[1]-w1[0];
	else if(telescopeParam==1)
		*g=max(w1[1]-w1[0], w1[degeneracy]-w1[degeneracy-1]);
	else
		*g=min(w1[1]-w1[0], w1[degeneracy]-w1[degeneracy-1]);

    free(R);
    free(psiPsiArrayN);
	free(psiHPsiArrayN);
	free(work);
	free(work2);
    free(psiHPsiArray);
    free(psiPsiArray);  

}	

double getError(double* errorMatrix, int dim){
	int n = dim, lda = n, info, lwork;
	double wkopt;
	double* work1;
	/* Local arrays */
	double w[n];
	/* Query and allocate the optimal workspace */
	lwork = -1;
	char mode;
	dsyev_("N", "Upper", &n, errorMatrix, &lda, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work1 = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_("N", "Upper", &n, errorMatrix, &lda, w, work1, &lwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues (getInternalGaps).\n" );
			exit( 1 );
	}	
	return fabs(w[0]);
}


/*
Used to find subspaces used for first excited states in 1st order tight binding/internal gaps as an error estimate in 0th order tight binding
*/

void getInternalGap(int wellWidth, double wellDepth, int bits, int wellNum, int subspace, double s, double* internalGaps, int* wellSubspaces, double* E0, double* E1){
    int dim=(bits+1-2*subspace);
	double* H=(double*)malloc(dim*dim*sizeof(double)); //Full H
	int x, y;
	//fill Operator Hamiltonian
    double sFac;
	if(wellNum==0){
		sFac=s;
	}else{
		sFac=1-s;
	}
	for(x=0; x<dim; x++){
		for(y=0; y<dim; y++){
			if(x==y){
				H[x*dim+y]=sFac*potOneWell(bits, wellWidth, wellDepth, x+subspace);			
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
	dsyev_("N", "Upper", &n, H, &lda, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work1 = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_("N", "Upper", &n, H, &lda, w, work1, &lwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues (getInternalGaps).\n" );
			exit( 1 );
	}	
	
    double gap;
    if(subspace==0){
		*E0=w[0];
		*E1=w[1];
        gap=w[1]-w[0]; 
        internalGaps[wellNum]=gap; 
        wellSubspaces[wellNum]=subspace;     
    }
    else{
		gap=w[0]-*E0; 
        if(gap<internalGaps[wellNum]){
			*E1=w[0];
            internalGaps[wellNum]=gap; 
            wellSubspaces[wellNum]=subspace;
        }  	
    }
    free(work1);
	free(H);	
}

/*
Used to generate array elements for tight binding
*/
//void generateArrayElements(double s,int bits,int wellNum1,int wellNum2,int wfNum1,int wfNum2,int* R,int numWells,int* wellWidths,double* wellDepths, int* wellSubspaces, double* psiPsi,double* psiHPsi, FILE* tbWfOutputFile, int getTbWfs){
void generateArrayElements(double s,int bits,int wellNum1,int wellNum2,int* R,int numWells, int* wellWidths,double* wellDepths,int* wellSubspaces,double** psiPsiArray,double** psiHPsiArray, FILE* tbWfOutputFile,int getTbWfs,int tbOrder, int numWfs, int* gammaArray, int* divisionIdx, double** errorMatrix){	
	int sep=R[wellNum1*numWells + wellNum2];
	int n1=sep;
	int n2=bits-sep;
	int dim1, dim2, subspace1, subspace2, dim;
	subspace1=wellSubspaces[wellNum1];	
	subspace2=wellSubspaces[wellNum2];
	int sig1, sig2;
	int i,j, well1Width, well2Width, h1, h2, x, y,k, gamma1, gamma2, gamma1p, gamma2p, h1p, h2p;
	double well1Depth, well2Depth, fac1, fac2;
	int Gamma1, Gamma2;
	double sum;
	int sigMax;
	if((subspace1==1||subspace2==1)&&tbOrder==1) sigMax=2;
	else sigMax=1;
	int u,v;

	for(u=0; u<sigMax; u++){
		for(v=0; v<sigMax; v++){
			if(u==1&&v==1) return; 
			if(subspace1!=subspace2&&(u!=0||v!=0)) return; //if subspaces are different and not in sig1=sig2=0 subspace relevant states won't have any overlap
			if(n1==0||n1==1){
				sig1=0; //this doesn't quite work fix this
				sig2=v;
				Gamma1=1; 
				Gamma2=(int) Gamma(n2, sig2);							
			}
			else if(n2==0||n2==1){
				sig1=u;
				sig2=0; 
				Gamma1=(int) Gamma(n1, sig1);
				Gamma2=1; 
				                                                                                                                                                   
			}else{
				sig1=u;
				sig2=v;
				Gamma1= (int) Gamma(n1, sig1);
				Gamma2= (int) Gamma(n2, sig2);				
			} 
			
			dim2=(n2+1-2*sig2)*Gamma2;
			dim1=(n1+1-2*sig1)*Gamma1;
			dim=dim1*dim2;
			
			well1Depth=wellDepths[wellNum1];
			well2Depth=wellDepths[wellNum2];
			well1Width=wellWidths[wellNum1];
			well2Width=wellWidths[wellNum2];
             double sFac1, sFac2;
			if(wellNum1==0){
				sFac1=(1-s);
			}else{
				sFac1=s;
			}
			if(wellNum2==0){
				sFac2=(1-s);
			}else{
				sFac2=s;
			}
			double* H1=(double *)malloc(dim*dim*sizeof(double)); //H just for 1st well
			double* H2=(double *)malloc(dim*dim*sizeof(double)); //H just for second well
			double* H=(double *)malloc(dim*dim*sizeof(double)); //H for full space
			double* HError=(double *)malloc(dim*dim*sizeof(double)); //H for full space - two wells (i.e. pertubation)
			//fill Operator Hamiltonian
			//printf("%i, %i, %i, %i, %i, %i\n",n1, n2, sig1, sig2, Gamma1, Gamma2);
			for(x=0; x<dim; x++){
				for(y=0; y<dim; y++){
					//INCLUDE GAMMA1 and GAMMA2 in here because raising is between same gamma 
					h1=x/(dim2*Gamma1)+sig1;
					gamma1=x/dim2%Gamma1;
					h2=x%(dim2)/Gamma2+sig2;
					gamma2=x%Gamma2;
					h1p=y/(dim2*Gamma1)+sig1;
					gamma1p=y/dim2%Gamma1;
					h2p=y%(dim2)/Gamma2+sig2;
					gamma2p=y%Gamma2;
					//printf("h1, gamma1, h2, gamma2: %i %i %i %i %i %i %i %i\n", h1, h2, gamma1, gamma2, h1p, h2p, gamma1p, gamma2p);
					if(h1==h1p&&h2==h2p&&gamma1==gamma1p&&gamma2==gamma2p){
						H[x*dim+y]=potAvg(h1, h2, n1, bits, wellWidths, wellDepths, wellNum1, wellNum2, numWells, R, s);	
						HError[x*dim+y]=potAvgErr(h1, h2, n1, bits, wellWidths, wellDepths, wellNum1, wellNum2, numWells, R,s);	
						H1[x*dim+y]=sFac1*pot(n1, n2, h1, h2, well1Width, 0.0, well1Depth, 0.0);
						H2[x*dim+y]=sFac2*pot(n1, n2, h1, h2, 0.0, well2Width, 0.0, well2Depth);
					}
					else if(h1==h1p&&h2+1==h2p&&gamma1==gamma1p&&gamma2==gamma2p&&n2!=0){
						H[x*dim+y]=(1-s)*hopUp(n2, h2, bits, sig2);
						HError[x*dim+y]=H[x*dim+y];
						H1[x*dim+y]=H[x*dim+y];
						H2[x*dim+y]=H[x*dim+y];
					}		
					else if(h1==h1p&&h2-1==h2p&&gamma1==gamma1p&&gamma2==gamma2p&&n2!=0){
						H[x*dim+y]=H[y*dim+x];
						HError[x*dim+y]=H[x*dim+y];
						H1[x*dim+y]=H[y*dim+x];
						H2[x*dim+y]=H[y*dim+x];
					}
					else if(h1+1==h1p&&h2==h2p&&gamma1==gamma1p&&gamma2==gamma2p&&n1!=0){
						H[x*dim+y]=(1-s)*hopUp(n1, h1, bits, sig1);
						HError[x*dim+y]=H[x*dim+y];
						H1[x*dim+y]=H[x*dim+y];
						H2[x*dim+y]=H[x*dim+y];
					}
					else if(h1-1==h1p&&h2==h2p&&gamma1==gamma1p&&gamma2==gamma2p&&n1!=0){
						H[x*dim+y]=H[y*dim+x];
						HError[x*dim+y]=H[y*dim+x];
						H1[x*dim+y]=H[y*dim+x];
						H2[x*dim+y]=H[y*dim+x];
					}
					else{
						H[x*dim+y]=0.0;
						HError[x*dim+y]=0.0;
						H1[x*dim+y]=0.0;
						H2[x*dim+y]=0.0;
					}			
				}		
			}

			// for(i=0; i<dim; i++){
			// 	for(j=0; j<dim; j++){
			// 		printf("% 2.3f ", H[i*dim+j]);
			// 	}
			// 	printf("\n");
			// }
			// printf("\n");

			int n = dim, lda = n, info, lwork;
			double wkopt;
			double* work1;
			double* work2;
			/* Local arrays */
			double w1[n];
			/* Query and allocate the optimal workspace */
			lwork = -1;
			dsyev_("V", "Upper", &n, H1, &lda, w1, &wkopt, &lwork, &info );
			lwork = (int)wkopt;
			work1 = (double*)malloc( lwork*sizeof(double) );
			/* Solve eigenproblem */
			dsyev_("V", "Upper", &n, H1, &lda, w1, work1, &lwork, &info );
			/* Check for convergence */
			if( info > 0 ) {
					printf( "The algorithm failed to compute eigenvalues (generateArrayElements(1)).\n" );
					exit( 1 );
			}
			
			double w2[n];
			/* Query and allocate the optimal workspace */
			lwork = -1;
			dsyev_("V", "Upper", &n, H2, &lda, w2, &wkopt, &lwork, &info );
			lwork = (int)wkopt;
			work2 = (double*)malloc( lwork*sizeof(double) );
			/* Solve eigenproblem */
			dsyev_("V", "Upper", &n, H2, &lda, w2, work2, &lwork, &info );
			/* Check for convergence */
			if( info > 0 ) {
					printf( "The algorithm failed to compute eigenvalues (generateArrayElements(2)).\n" );
					exit( 1 );
			}
			//Now get all the matrix elements we need
			//start with GS and 1ex in 0,0 space
			int count,init;
			if(sig1==0&&sig2==0){
				init=0;			
				if(tbOrder==0)
					count=1;
				else	
					count=2; //gs plus first excited state in 0,0 subspace
			}else{
				if(sig1==1||sig2==1){
					init=2;
					count=gammaArray[wellNum1]-1; //should equal gammaArray[wellNum2] or else we already returned
					//gammaArray gives number of excited states plus 1 gs?
				}
				else{
					init=0;
					count=0;//should never get here
				}
			} 
			double psi1[dim];
			double psi2[dim];
			double Hpsi2[dim];
			double HErrorpsi2[dim];
			for(x=0; x<count; x++){
				for(y=0; y<count; y++){
					fac1=1;
					fac2=1;
					if(H1[dim*x]<2e-16) fac1=-1;
					if(H2[dim*y+n1*dim2*Gamma1]<2e-16) fac2=-1;
					for(i=0; i<dim; i++){
						psi1[i]=fac1*H1[dim*x+i];
						psi2[i]=fac2*H2[dim*y+i];
					}
					if(getTbWfs==1){
						//write out wfs
						fprintf(tbWfOutputFile, "%1.15f, %i, %i, %i, %i, %i, %i, %i, %i, %i, %i\n", s, bits, n1, n2,dim, dim1, dim2, wellNum1, x, wellNum2, y);
						for(i=0; i<dim; i++){
							fprintf(tbWfOutputFile, "% 1.15f, ", psi1[i]);
						}
						fprintf(tbWfOutputFile, "\n");
						for(i=0; i<dim; i++){
							fprintf(tbWfOutputFile, "% 1.15f, ", psi2[i]);
						}
						fprintf(tbWfOutputFile, "\n");
					}
					
					//Get H times state 1						
					for(i=0; i<dim; i++){
						sum=0.0;
						for(j=0; j<dim; j++){
							sum+=H[i*dim+j]*psi2[j];			
						}
						Hpsi2[i]=sum;		
					}
						
					//get H matrix element 
					sum=0.0;						
					for(i=0; i<dim; i++){
						sum+=psi1[i]*Hpsi2[i];					
					}
					

					(*psiHPsiArray)[(divisionIdx[wellNum1]+(x+init))*numWfs+(divisionIdx[wellNum2]+(y+init))]=sum;
					(*psiHPsiArray)[(divisionIdx[wellNum2]+(y+init))*numWfs+(divisionIdx[wellNum1]+(x+init))]=sum;

					//Get HError times state 1						
					for(i=0; i<dim; i++){
						sum=0.0;
						for(j=0; j<dim; j++){
							sum+=HError[i*dim+j]*psi2[j];			
						}
						HErrorpsi2[i]=sum;		
					}
						
					//get H matrix element 
					sum=0.0;						
					for(i=0; i<dim; i++){
						sum+=psi1[i]*HErrorpsi2[i];					
					}
					

					(*errorMatrix)[(divisionIdx[wellNum1]+(x+init))*numWfs+(divisionIdx[wellNum2]+(y+init))]=sum;
					(*errorMatrix)[(divisionIdx[wellNum2]+(y+init))*numWfs+(divisionIdx[wellNum1]+(x+init))]=sum;


					//get S matrix element 
					if(wellNum1==wellNum2){
						
						if(x==y){
							(*psiPsiArray)[(divisionIdx[wellNum1]+(x+init))*numWfs+(divisionIdx[wellNum2]+(y+init))]=1.0;
						}else{
							(*psiPsiArray)[(divisionIdx[wellNum1]+(x+init))*numWfs+(divisionIdx[wellNum2]+(y+init))]=0.0;
							(*psiPsiArray)[(divisionIdx[wellNum2]+(y+init))*numWfs+(divisionIdx[wellNum1]+(x+init))]=0.0;
						}
					}else{
						sum=0.0;
						for(i=0; i<dim; i++){
							sum+=psi1[i]*psi2[i];					
						}
						if(fabs(sum)<2e-16) sum=0.0;
						(*psiPsiArray)[(divisionIdx[wellNum1]+(x+init))*numWfs+(divisionIdx[wellNum2]+(y+init))]=sum;
						(*psiPsiArray)[(divisionIdx[wellNum2]+(y+init))*numWfs+(divisionIdx[wellNum1]+(x+init))]=sum;
					}

				}
			}
			free(H1);
			free(H2);
			free(H);
			free( (void*)work1 ); 
			free( (void*)work2 );	
		}
	}
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
Get potential for one well of fixed width
*/
double potOneWell(int bits,int width,double depth, int w){
	if(w<width) return depth;
	return 0.0;
}


/*
Potential term in 2 well basis
*/
double potAvg(int h1,int h2, int n1, int bits, int* wellWidths,double* wellDepths, int i, int j, int numWells, int* R, double s){
	//note must guarantee wells don't overlap beforehand
	int x;
	double Vtot=0.0;
	//printf("i,j: %i, %i\n", i, j);
	if(i==0 && j==0) numWells=2; //prevents from multiply counting potential in case of a non-well which is always placed first
	for(x=0; x<numWells; x++){
		if(x==i){
			if(h1+h2<wellWidths[i]){
				if(x==0){
					return s*wellDepths[i];
				}else{					
					return (1-s)*wellDepths[i];
				}
			}
		}
		else if(x==j){
			if(h2+(n1-h1)<wellWidths[j]){
				if(x==0){
					return s*wellDepths[j];
				}else{					
					return (1-s)*wellDepths[j];
				}
			}
		}	
		else{			
			if(x==0){
				Vtot+=N(R[i*numWells+x], R[j*numWells+x], n1, h1, h2, wellWidths[x]-1, bits)*wellDepths[x]*s; 
			}else{					
				Vtot+=N(R[i*numWells+x], R[j*numWells+x], n1, h1, h2, wellWidths[x]-1, bits)*wellDepths[x]*(1-s); 
			}	
		//	printf("%i\n", R[i*numWells+x]);			
				//N gets the overlap of the shell w/ radius h1+h2 and the well x			
		}
	}
	return Vtot/sqrt(binomial(n1, h1)*binomial(bits-n1, h2));	
}


/*
Potential term in 2 well basis
*/
double potAvgErr(int h1,int h2, int n1, int bits, int* wellWidths,double* wellDepths, int i, int j, int numWells, int* R, double s){
	//note must guarantee wells don't overlap beforehand
	int x;
	double Vtot=0.0;
	//printf("i,j: %i, %i\n", i, j);
	if(i==0 && j==0) numWells=2; //prevents from multiply counting potential in case of a non-well which is always placed first
	for(x=0; x<numWells; x++){
		if(x==i){
			if(h1+h2<wellWidths[i]){
				return 0.0;
			}
		}
		else if(x==j){
			if(h2+(n1-h1)<wellWidths[j]){
				return 0.0;
			}
		}	
		else{			
			if(x==0){
				Vtot+=N(R[i*numWells+x], R[j*numWells+x], n1, h1, h2, wellWidths[x]-1, bits)*wellDepths[x]*s; 
			}else{					
				Vtot+=N(R[i*numWells+x], R[j*numWells+x], n1, h1, h2, wellWidths[x]-1, bits)*wellDepths[x]*(1-s); 
			}	
		//	printf("%i\n", R[i*numWells+x]);			
				//N gets the overlap of the shell w/ radius h1+h2 and the well x			
		}
	}
	return Vtot/sqrt(binomial(n1, h1)*binomial(bits-n1, h2));
}


/*Returns number of points of overlap of 2 hamming spheres with one hamming ball*/
int N(int R12, int R23, int n1, int h1, int h2, int w, int bits){
	int r1=h1+h2;
	int r2=h2+(n1-h1);
	int r3=0;
	int x=0.5*(n1+R12-R23);
	int y=0.5*(n1-R12+R23);
	int z=0.5*(-n1+R12+R23);
	int overlaps=0;
	double kx;
	double ky, kz, kl;
	for(kx=0; kx<=x; kx++){
		for(r3=0; r3<=w; r3++){
			ky=0.5*(y+r1-r2-x)+kx;
			kz=0.5*(z+r1-r3-x)+kx;
			kl=-0.5*(y+z-r2-r3)-kx;
			if(ky>=0&&ky<=y&&kz>=0&&kz<=z&&kl>=0&&kl<=(bits-x-y-z)&&ceilf(ky) == ky&&ceilf(kz) == kz&&ceilf(kl) == kl) overlaps+=binomial(x,kx)*binomial(y,ky)*binomial(z,kz)*binomial(bits-x-y-z,kl);
		}		
	}
	//printf("N, h1, h2: %i, %i, %i\n", overlaps, h1, h2);
	return overlaps;
}

/*returns the binomial(n,k)*/
double binomial(int n, int k) {
	if(n==0&&k==0) return 1;
	else if(k<0) return 0;
	else if(n<0) return 1;
	else if(n==0&&k!=0) return 0;

  double numerator;
  double denominator;
  int i;
  if(k>n) return 0;
  if(k > n/2) k = n-k;
  numerator = 1.0;
  denominator = 1.0;
  for(i = n-k+1; i <= n; i++) numerator *= (double)i;
  for(i = 2; i <= k; i++) denominator *= (double)i;
  return numerator/denominator;
}

double Gamma(int bits,int sig){
	return binomial(bits, sig)-binomial(bits, sig-1);
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
