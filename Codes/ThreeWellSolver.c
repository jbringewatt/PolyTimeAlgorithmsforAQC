/***********************************
************************************
Author: Jake Bringewatt

Use: Exact three well solver. Only gets energies in sig_1, sig_2, sig_3, sig_4=0 subspace.

Dependencies: Lapack, Blas

Compile: gcc -o ThreeWellSolver ./ThreeWellSolver-lm -lblas -llapack
Run: ./ThreeWellSolver wellInputFile paramInputFile

Input: 
Standard well input file (see README) [wellInputFile]:
Has the standard form for all codes:
number of qubits
number of wells (will return with error if not 3)
well depth (double)(<0)
well width (integer)(in [1,n])

e.g.
10
3
1111011101
-2.1
1
0000101110
-2
2
1110101110
-2
1

One well solver parameter input file [paramInputFile]:
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

Output: 2 csv files:
1. runname.tw.energies.output
Gives the energies and gaps at each s where each row is
s, E0 (sigma=0), E1 (sigma=0), E1*(sigma=1), E1-E0, E1*-E0

2. runname.tw.wfs.output
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
void getStatesHelper(char** wellLocations, int* wellWidths, double* wellDepths, int bits, int numWells,double s, double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs);
extern void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info);
double pot(int h1,int h2,int h3, int h4, int n1, int n2, int n3, int well1Width,int well2Width,int well3Width, double well1Depth,double well2Depth, double well3Depth, double s);
double hopUp(int n,int h,int bits);
double hopDown(int n,int h,int bits);
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
  
  if(numWells!=3){
      printf("ERROR: Not three wells in input file.\n");
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
		if(wellLocations[0][i]=='1'&&wellLocations[1][i]=='1'&&wellLocations[2][i]=='1'){
			wellLocations[0][i]='0';
			wellLocations[1][i]='0';
			wellLocations[2][i]='0';
		}
		else if(wellLocations[0][i]=='0'&&wellLocations[1][i]=='1'&&wellLocations[2][i]=='1'){
			wellLocations[0][i]='1';
			wellLocations[1][i]='0';
			wellLocations[2][i]='0';
		}
		else if(wellLocations[0][i]=='1'&&wellLocations[1][i]=='0'&&wellLocations[2][i]=='1'){
			wellLocations[0][i]='0';
			wellLocations[1][i]='1';
			wellLocations[2][i]='0';
		}
		else if(wellLocations[0][i]=='1'&&wellLocations[1][i]=='1'&&wellLocations[2][i]=='0'){
			wellLocations[0][i]='0';
			wellLocations[1][i]='0';
			wellLocations[2][i]='1';
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
void getStates(double s,int bits, int numWells, char** wellLocations,int* wellWidths,double* wellDepths, double* g, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs, int telescopeParam, double*E0){
	printf("s: %f\n", s);
   
    double gaps[4];	
	//check only sig1=0, sig2=0, sig3=0, sig4=0 subspace - this is the only relevant gap for three wells, simplifies the code
	getStatesHelper(wellLocations, wellWidths, wellDepths, bits, numWells, s, gaps, E0, energyOutputFile, wfOutputFile, getWfs);
     
    	
	if(telescopeParam==0)
		*g=gaps[0];
	else if(telescopeParam==1)
		*g=max(gaps[0], gaps[0]);
	else
		*g=min(gaps[0], gaps[0]);		
}	


/*
Used to find energies and wfs for a given subspace
*/
void getStatesHelper(char** wellLocations, int* wellWidths, double* wellDepths, int bits, int numWells,double s, double* gaps, double* E0, FILE* energyOutputFile, FILE* wfOutputFile, int getWfs){
	//find well separations
	int i, j, k, n1, n2, n3, n4;
	//Find n1
	n1=0;
	n2=0;	
	n3=0;
	n4=0;
	for(k=0; k<bits; k++){		
		if(wellLocations[0][k]=='1') n1++;
		else if(wellLocations[1][k]=='1') n2++;
		else if(wellLocations[2][k]=='1') n3++;	
		else n4++;
	}
//	printf("n1,n2,n3, n4: %i, %i, %i, %i\n", n1,n2,n3, n4);
	int dim=(n1+1)*(n2+1)*(n3+1)*(n4+1);
	double well1Depth=wellDepths[0];
	double well2Depth=wellDepths[1];
	double well3Depth=wellDepths[2];
	int well1Width=wellWidths[0];
	int well2Width=wellWidths[1];
	int well3Width=wellWidths[2];
	double* H=(double *)malloc(dim*dim*sizeof(double));
	int h1=0, h2=0, h3=0, h4=0;
	int x, y;
	//fill Operator Hamiltonian
	for(x=0; x<dim; x++){
		for(y=0; y<dim; y++){
			h1=x/((n2+1)*(n3+1)*(n4+1));
			h2=x/((n3+1)*(n4+1))%(n2+1);
			h3=x/(n4+1)%(n3+1);
			h4=x%(n4+1);
			
			//printf("h1, h2, h3, h4: %i, %i, %i, %i\n", h1, h2, h3, h4);
			
			if(x==y){
				H[x*dim+y]=pot(h1, h2, h3, h4, n1, n2, n3, well1Width, well2Width, well3Width, well1Depth, well2Depth, well3Depth,s);					
			}
			else if(y==x+1&&n4!=0){
				H[x*dim+y]=(1-s)*hopUp(n4, h4, bits);				
			}		
			else if(y==x-1&&n4!=0){
				H[x*dim+y]=(1-s)*hopDown(n4, h4, bits);				
			}
			else if(y==x+n4+1&&n3!=0){
				H[x*dim+y]=(1-s)*hopUp(n3, h3, bits);				
			}
			else if(y==x-n4-1&&n3!=0){
				H[x*dim+y]=(1-s)*hopDown(n3, h3, bits);
			}
			else if(y==x+((n3+1)*(n4+1))&&n2!=0){				
				H[x*dim+y]=(1-s)*hopUp(n2, h2, bits);								
			}
			else if(y==x-((n3+1)*(n4+1))&&n2!=0){
				H[x*dim+y]=(1-s)*hopDown(n2, h2, bits);				
			}		
			else if(y==x+((n2+1)*(n3+1)*(n4+1))&&n1!=0){	
				H[x*dim+y]=(1-s)*hopUp(n1, h1, bits);								
			}
			else if(y==x-((n2+1)*(n3+1)*(n4+1))&&n1!=0){
				H[x*dim+y]=(1-s)*hopDown(n1, h1, bits);				
			}					
			else{
				H[x*dim+y]=0.0;				
			}			
		}		
	}
	
	
/* 	 printf("H\n");
	 for(x=0; x<dim; x++){
	  for(y=0; y<dim; y++){
		  printf("%f ", H[x*dim+y]);		  
		}
		printf("\n");
	}    */ 
	
	//once you have TBH for a given s, calculate the eigenvalue gap
	/* Locals */
	int n = dim, lda = n, info, lwork;
	double wkopt;
	double* work;
	/* Local arrays */
	double w[n];
	lwork = -1;
	dsyev_( "Vectors", "Upper", &n, H, &lda, w, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Solve eigenproblem */
	dsyev_( "Vectors", "Upper", &n, H, &lda, w, work, &lwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
			printf( "The algorithm failed to compute eigenvalues.\n" );
			exit( 1 );
	}
  
		*E0=w[0];
    fprintf(energyOutputFile, "%1.15f, %1.15f, %1.15f, ", s, w[0], w[1]);
		gaps[0]=w[1]-w[0];       
    fprintf(energyOutputFile, "%1.15f, %1.15f\n", w[0], gaps[0]);	
    
    
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
    	if(H[dim]>=-1e-15) fac1=1.0;
			else fac1=-1.0;		
			for(j=0; j<n; j++){
				fprintf(wfOutputFile, "% 1.15f, ", fac1*H[dim+j]);
			}
			fprintf(wfOutputFile, "\n"); 
  }          
	free(work);
	free(H);	
}

/*
3 well potential
*/
double pot(int h1,int h2,int h3, int h4, int n1, int n2, int n3, int well1Width,int well2Width,int well3Width, double well1Depth,double well2Depth, double well3Depth, double s){
	//note must guarantee wells don't overlap beforehand
	if(h2+h3+h4+(n1-h1)<well1Width){
		return well1Depth*s;
	}
	else if(h1+(n2-h2)+h3+h4<well2Width){
		return well2Depth*s;
	}
	else if(h1+h2+(n3-h3)+h4<well3Width){
		return well3Depth*s;
	}
	return 0;	
}

double hopUp(int n,int h,int bits){
	double dn, dh, dbits;
	dn = (double) n;
	dh = (double) h;
	dbits = (double) bits;
	return -1*sqrt((dh+1.0)*(dn-dh))/dbits;
}

double hopDown(int n,int h,int bits){
	double dn, dh, dbits;
	dn = (double) n;
	dh = (double) h;
	dbits = (double) bits;
	return -1*sqrt(dh*(dn-dh+1.0))/dbits;
}