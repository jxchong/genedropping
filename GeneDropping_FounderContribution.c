#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MersenneTwister.h"

MTRand twist;

// integer from 0 to n-1
int randInt(int max) { 
    //return (int)((float) max*rand()/(RAND_MAX+1.0));
	return twist.randInt(max-1);
}

// float in range
float randRange(int a, float b) {
	//return ( (b-a) * ((float) rand() / (RAND_MAX+1.0)) )+a;
	return twist.rand(b-a) + a;
}

// do the Allele selection
int randAllele(int numMom, int numDad, float selection) {
	if ((numMom == 0) && (numDad == 0)) {
		return 0;
	}
	else if ((numMom == 1) && (numDad == 0)) {
		return randInt(2);
	}
	else if ((numMom == 2) && (numDad == 0)) {
		return 1;
	}
	else if ((numMom == 0) && (numDad == 1)) {
		return randInt(2);
	}
	else if ((numMom == 1) && (numDad == 1)) {
		float val = randRange(0,(3+selection));
		if (val < 1.0) {
			return 0;
		}
		else if (val < 3.0) {
			return 1;
		}
		else {
			return 2;
		}
	}
	else if ((numMom == 2) && (numDad == 1)) {
		float val = randRange(0,(1+selection));
		if (val < 1.0) {
			return 1;
		}
		else {
			return 2;
		}
	}
	else if ((numMom == 0) && (numDad == 2)) {
		return 1;
	}
	else if ((numMom == 1) && (numDad == 2)) {
		float val = randRange(0,(1+selection));
		if (val < 1.0) {
			return 1;
		}
		else {
			return 2;
		}
	}
	else if ((numMom == 2) && (numDad == 2)) {
		return 2;
	}
}


struct genocount {
	int AA, Aa, aa;
};

void init_genocount (struct genocount* thisinstance) {
	thisinstance->AA = 0;
	thisinstance->Aa = 0;
	thisinstance->aa = 0;
}

int main(int argc, char** argv) {
	// test for correct inputs
	if (argc != 2) {
		printf("Usage: ./GeneDropping_CarrierBurden [config file]\n\n");
		printf("Config file format (one per line): \n 1) path/name of pedigree file\n 2) number of hutterites in pedigree\n 3) number of founders in pedigree\n 4) path/name of cohort file\n 5) number of cohorts in the file (assumes that 0 counts as a group containing all no-cohort findivs)\n 6) number of sims desired per founder\n 7) fitness penalty of mutations (0 to 1)\n 8) number of mutations per founder 9) path/name of output file\n\n");
		exit(1);
	}
	
	// load information from config file
	int numHut, numFounder, num_cohorts, num_desired, selection, numFounderMutations;
    
	// char *pedfile; 							// initialize pointers to files
	// char *cohortfile; 						// initialize pointers to files
	// char *outputfile; 						// initialize pointers to files
	char pedfile[2000];
	char cohortfile[2000];
	char outputfile[2000];
	// printf("Loading information from config file\n");			// DEBUG
	
	FILE* configfile = fopen(argv[1], "r");
	if (configfile == NULL) {
		printf("ERROR: Config file does not exist\n");
		exit(1);
	}
	fscanf(configfile, "%s", pedfile);
	printf("Pedigree file = %s\n", pedfile);			// DEBUG
	fscanf(configfile, "%i", &numHut);
	printf("Number of Hutterites in pedigree = %i\n", numHut);			// DEBUG
	fscanf(configfile, "%i", &numFounder);
	printf("Number of founders in pedigree = %i\n", numFounder);			// DEBUG
	fscanf(configfile, "%s", cohortfile);
	printf("Cohort file = %s\n", cohortfile);			// DEBUG
	fscanf(configfile, "%i", &num_cohorts);
	printf("Number of cohorts = %i\n", num_cohorts);			// DEBUG
	fscanf(configfile, "%i", &num_desired);
	printf("Number of simulations per founder = %i\n", num_desired);			// DEBUG
	fscanf(configfile, "%i", &selection);
	printf("Selection/fitness value = %i\n", selection);			// DEBUG
    fscanf(configfile, "%i", &numFounderMutations);
    printf("Number of founder mutations = %i\n", numFounderMutations);
	fscanf(configfile, "%s", outputfile);
	printf("Output file = %s\n\n", outputfile);			// DEBUG
	fclose(configfile);
	
	// load the pedigree
	FILE* pedfile_handle = fopen(pedfile, "r");
	if (pedfile_handle == NULL) {
		printf("ERROR: No Such Pedigree\n");
		exit(1);
	}
	
	// printf("Loading pedigree file.\n");					// DEBUG
	int findivs[numHut];
	int numAllele[numHut];
	int* toMom[numHut];
	int* toDad[numHut];
    
    int binHeights[numFounder*2*numFounderMutations];
    int carrierBurden[numHut];                               // track number of alleles per person, over a single iteration (from 1:64 founders)
    int numMutationsSurvive[numFounder*numFounderMutations];
    float meannumMutationsSurvive = 0;
    int nMutPerFounderSurvive[numFounder];                 // track number of mutations per founder surviving in current population
    float meanPropSurvivingMutPerFounder[numFounder];           // expected proportion of surviving mutations inherited from a single founder
    	
	int i, j;
	for (i = 0; i < numHut; i++) {
		int findiv, m_findiv, f_findiv;
		fscanf(pedfile_handle, "%i\t%i\t%i", &findiv, &m_findiv, &f_findiv);
		findivs[i] = findiv;
		toMom[i] = NULL;
		toDad[i] = NULL;
		for (j = 0; j < i; j++ ) {
			if((m_findiv != 0) && (findivs[j] == m_findiv)) {
				toMom[i] = &(numAllele[j]);
			}
			if((f_findiv != 0) && (findivs[j] == f_findiv)) {
				toDad[i] = &(numAllele[j]);
			}
		}
	}
	
	fclose(pedfile_handle);
	
	
	// load the cohorts
	FILE* cohortfile_handle = fopen(cohortfile, "r");
	if (cohortfile_handle == NULL) {
		printf("ERROR: Birth cohort file doesn't exist\n");
		exit(1);
	}
	// printf("Loading cohorts.\n");						// DEBUG

	struct genocount cohort_counts[num_cohorts];			// store genotype counts for each cohort, where each genocount is an element in the array cohort_counts

	// printf("Initialize cohort_counts.\n");						// DEBUG	
	for (i = 0; i<num_cohorts; i++) {
		init_genocount( &cohort_counts[i] );
	}
	
	// printf("Zeroing out to_cohort.\n");						// DEBUG
	struct genocount* to_cohort[numHut];
	for (i=0; i<numHut; i++) {						// zero out pointers
		to_cohort[i] = NULL;
	}
	
	int myfindiv, mycohort;
	int totalsubjCounted = 0;
	while (fscanf(cohortfile_handle, "%d\t%d", &myfindiv, &mycohort) != EOF) {
		if (mycohort != 0) {
		    totalsubjCounted++;
			for (i = 0; i < numHut; i++) {
				if (findivs[i] == myfindiv) {
					to_cohort[i] = &( cohort_counts[mycohort] );		// create pointer to
				}
			}
		}
	}
	fclose(cohortfile_handle);

	
	// FINISHED LOADING DATA AND SETTING UP STRUCTURES
	// BEGIN SIMULATION
	
	printf("\nBegin simulation\n");
	
	FILE* outfile_handle = fopen(outputfile, "w");
	if (outfile_handle == NULL) {
		printf("ERROR: Could Not Open Output File\n");
		exit(1);
	}
	
	// printf("set num_trials val.\n");				// DEBUG
	int num_trials = 0;
	for (num_trials=0; num_trials<num_desired; num_trials++) {
		printf("Trial #%d\n", num_trials);
		
		// initialize total carrier burden and count of surviving founder mutations
		for (i = 0; i < numHut; i++) {
			carrierBurden[i] = 0;
		}
	    for (i = 0; i < numFounder*numFounderMutations; i++) {
			numMutationsSurvive[i] = 0;
		}
	    for (i = 0; i < numFounder; i++) {
			nMutPerFounderSurvive[i] = 0;
		}
	
    	int chosenFounderidx = 0;
    	for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {			// iterate over all founders
    		int chosenFounder = findivs[chosenFounderidx];
	        
            int nmut = 0;
            for (nmut=0; nmut<numFounderMutations; nmut++) {
    			// initialize and set founder alleles
    			for (i = 0; i < numFounder; i++) {
    				numAllele[i] = 0;
    			}
    			numAllele[chosenFounderidx] = 1;		
		
    			// simulate the rest
    			for (i = numFounder; i < numHut; i++) {
    				int numMom = *(toMom[i]);
    				int numDad = *(toDad[i]);			
    				numAllele[i] = randAllele(numMom, numDad, selection);
    			}
    			
    			// finished dropping alleles for this founder, save information
    			for (j = 0; j < numHut; j++) {
                    carrierBurden[j] += numAllele[j];
                    
                    // determine if founder mutation "survives" to desired cohort
                    if (to_cohort[j] == &(cohort_counts[1])) {
                        numMutationsSurvive[chosenFounderidx+numFounder*nmut] += numAllele[j];
                    }
    			}
			}
		}        
		// finish of a single iteration over all founders, count number of subjects falling into a given carrier burden bin
		for (i = 0; i < numHut; i++) {
		    if (to_cohort[i] == &(cohort_counts[1])) {
                binHeights[carrierBurden[i]]++;
            }
		}
		    
        int nmut = 0;        
        // for (i = 0; i < numFounder*numFounderMutations; i++) {        
        for (nmut=0; nmut<numFounderMutations; nmut++) {    
        	for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {
                if (numMutationsSurvive[chosenFounderidx+numFounder*nmut] > 0) {
                    meannumMutationsSurvive++;
                    nMutPerFounderSurvive[chosenFounderidx]++;
                }
        	}
    	}
    	
    	int totMutSurvive = 0;
    	for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {
            totMutSurvive += nMutPerFounderSurvive[chosenFounderidx];
        }
        
        for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {
           meanPropSurvivingMutPerFounder[chosenFounderidx] += (1.0*nMutPerFounderSurvive[chosenFounderidx])/totMutSurvive;
        }
        
    	
	}  // finish all iterations/trials




	// calculate mean (over all iterations bin height)
    float meanbinHeights[numFounder*2*numFounderMutations];
    for (i = 0; i < numFounder*2*numFounderMutations; i++) {
        meanbinHeights[i] = (1.0*binHeights[i])/num_trials;
    }
    
    // calculate mean number of surviving founder mutations
    meannumMutationsSurvive = meannumMutationsSurvive/num_trials;
    
    // calculate mean proportion of surviving mutations attributable to each founder
    int chosenFounderidx = 0;
    for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {
         meanPropSurvivingMutPerFounder[chosenFounderidx] = meanPropSurvivingMutPerFounder[chosenFounderidx]/num_trials;
     }
	
	// print results to file							// DEBUG
	// printf("print results to output\n");					// DEBUG
	
	fprintf(outfile_handle, "n mutations per founder = %i\n", numFounderMutations);	
	fprintf(outfile_handle, "fitness of homozygotes = %i\n", selection);	
    fprintf(outfile_handle, "mean surviving founder mutations = %f\n", meannumMutationsSurvive);			
	fprintf(outfile_handle, "initial total founder mutations = %i\n\n", numFounder*numFounderMutations);
	
    fprintf(outfile_handle, "proportion of surviving mutations from founder x\n");
    
    for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {
        fprintf(outfile_handle, "%i\t%f\n", findivs[chosenFounderidx], meanPropSurvivingMutPerFounder[chosenFounderidx]);	
    }
	
	
    //     fprintf(outfile_handle, "burden\tnsubj\tpropsubj\n");            // DEBUG
    // 
    // for (i = 0; i < numFounder*2*numFounderMutations; i++) {
    //     if (meanbinHeights[i] > 0.00000001) {
    //             float propsubj = meanbinHeights[i]/totalsubjCounted; 
    //         fprintf(outfile_handle, "%d\t%f\t%f\n", i, meanbinHeights[i], propsubj);
    //     }
    //     }    
    // 
	fclose(outfile_handle);	
	exit(0);
}

