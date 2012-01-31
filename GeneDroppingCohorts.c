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
		printf("Usage: ./GeneDroppingCohorts [config file]\n\n");
		printf("Config file format (one per line): \n 1) path/name of pedigree file\n 2) number of hutterites in pedigree\n 3) number of founders in pedigree\n 4) path/name of cohort file\n 5) number of cohorts in the file (assumes that 0 counts as a group containing all no-cohort findivs)\n 6) number of sims desired per founder\n 7) fitness penalty of mutations (0 to 1)\n 8) path/name of output file\n\n");
		exit(1);
	}
	
	// load information from config file
	int numHut, numFounder, num_cohorts, num_desired, selection;
	// char *pedfile; 							// initialize pointers to files
	// char *cohortfile; 						// initialize pointers to files
	// char *outputfile; 						// initialize pointers to files
	char pedfile[2000];
	char cohortfile[2000];
	char outputfile[2000];
	char outputfileprophet[2000];
	char outputfileprophom[2000];
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
	fscanf(configfile, "%s", outputfile);
	printf("Output file = %s\n", outputfile);			// DEBUG
	fscanf(configfile, "%s", outputfileprophet);
	printf("Output file = %s\n", outputfileprophet);			// DEBUG
	fscanf(configfile, "%s", outputfileprophom);
	printf("Output file = %s\n\n", outputfileprophom);			// DEBUG
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
	for (i = 0; i < num_cohorts; i++) {
		init_genocount( &cohort_counts[i] );
	}
	
	// printf("Zeroing out to_cohort.\n");						// DEBUG
	struct genocount* to_cohort[numHut];
	for (i=0; i<numHut; i++) {						// zero out pointers
		to_cohort[i] = NULL;
	}
	
	int myfindiv, mycohort;
	while (fscanf(cohortfile_handle, "%d\t%d", &myfindiv, &mycohort) != EOF) {
		if (mycohort != 0) {
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
	FILE* outfileprophom_handle = fopen(outputfileprophom, "w");
	if (outfileprophom_handle == NULL) {
		printf("ERROR: Could Not Open Output File\n");
		exit(1);
	}
	FILE* outfileprophet_handle = fopen(outputfileprophet, "w");
	if (outfileprophet_handle == NULL) {
		printf("ERROR: Could Not Open Output File\n");
		exit(1);
	}
	
	fprintf(outfile_handle, "founder\tcohort1freq\tcohort2freq\tcohort3freq\tcohort4freq\tcohort5freq\n");  			// DEBUG
	fprintf(outfileprophet_handle, "founder\tcohort1freq\tcohort2freq\tcohort3freq\tcohort4freq\tcohort5freq\n");  			// DEBUG
	fprintf(outfileprophom_handle, "founder\tcohort1prophom\tcohort2prophom\tcohort3prophom\tcohort4prophom\tcohort5prophom\n");
	
	int chosenFounderidx = 0;
	for (chosenFounderidx=0; chosenFounderidx<numFounder; chosenFounderidx++) {			// iterate over all founders
		int chosenFounder = findivs[chosenFounderidx];
		printf("Simulating for founder #%d %d\n", chosenFounderidx, chosenFounder);
	
	
		// printf("set num_trials val.\n");				// DEBUG
		int num_trials = 0;
		for (num_trials=0; num_trials<num_desired; num_trials++) {
			//set founder alleles
			// printf("set founders.\n");				// DEBUG
			for (i = 0; i < numFounder; i++) {
				numAllele[i] = 0;
			}
			numAllele[chosenFounderidx] = 1;
		
			// printf("set others.\n");				// DEBUG
		
			// simulate the rest
			for (i = numFounder; i < numHut; i++) {
				/*if ((toMom[i]) == NULL) {
					printf("no mom! %d\n", findivs[i]);
					exit(1);
				}
				if ((toDad[i]) == NULL) {
					printf("no dad! %d\n", findivs[i]);
					exit(1);
				}*/
				int numMom = *(toMom[i]);
				//printf("mom's stuff : %d\n", numMom);
				int numDad = *(toDad[i]);
				//printf("dad's stuff : %d\n", numMom);
			
				numAllele[i] = randAllele(numMom, numDad, selection);
			}
		
			// printf("storing data.\n");				\\ DEBUG
		
			// finished dropping alleles for this iterations, save informatino
			// calculate frequency in current cohort
			for (j = 0; j < numHut; j++) {
				int findiv = findivs[j];
				int mycount = numAllele[j];
				if (to_cohort[j] != NULL) {
					if (mycount == 0) {
						to_cohort[j]->AA += 1;
					}
					else if (mycount == 1) {
						to_cohort[j]->Aa += 1;
					}
					else if (mycount == 2) {
						to_cohort[j]->aa += 1;
					}
				}
			}
		
			// print results to file							// DEBUG
			// printf("print results to output\n");					// DEBUG
			float freq1 = (cohort_counts[1].Aa + 2.0 * cohort_counts[1].aa)/(2.0*(cohort_counts[1].AA + cohort_counts[1].Aa + cohort_counts[1].aa));
			float freq2 = (cohort_counts[2].Aa + 2.0 * cohort_counts[2].aa)/(2.0*(cohort_counts[2].AA + cohort_counts[2].Aa + cohort_counts[2].aa));
			float freq3 = (cohort_counts[3].Aa + 2.0 * cohort_counts[3].aa)/(2.0*(cohort_counts[3].AA + cohort_counts[3].Aa + cohort_counts[3].aa));
			float freq4 = (cohort_counts[4].Aa + 2.0 * cohort_counts[4].aa)/(2.0*(cohort_counts[4].AA + cohort_counts[4].Aa + cohort_counts[4].aa));
			float freq5 = (cohort_counts[5].Aa + 2.0 * cohort_counts[5].aa)/(2.0*(cohort_counts[5].AA + cohort_counts[5].Aa + cohort_counts[5].aa));
			fprintf(outfile_handle, "%d\t%f\t%f\t%f\t%f\t%f\n", chosenFounder, freq1, freq2, freq3, freq4, freq5); 

			float hetfreq1 = (1.0 * cohort_counts[1].Aa)/((cohort_counts[1].AA + cohort_counts[1].Aa + cohort_counts[1].aa));
			float hetfreq2 = (1.0 * cohort_counts[2].Aa)/((cohort_counts[2].AA + cohort_counts[2].Aa + cohort_counts[2].aa));
			float hetfreq3 = (1.0 * cohort_counts[3].Aa)/((cohort_counts[3].AA + cohort_counts[3].Aa + cohort_counts[3].aa));
			float hetfreq4 = (1.0 * cohort_counts[4].Aa)/((cohort_counts[4].AA + cohort_counts[4].Aa + cohort_counts[4].aa));
			float hetfreq5 = (1.0 * cohort_counts[5].Aa)/((cohort_counts[5].AA + cohort_counts[5].Aa + cohort_counts[5].aa));
			fprintf(outfileprophet_handle, "%d\t%f\t%f\t%f\t%f\t%f\n", chosenFounder, freq1, freq2, freq3, freq4, freq5);

            float prophom1 = (1.0 * cohort_counts[1].aa/(cohort_counts[1].AA + cohort_counts[1].Aa + cohort_counts[1].aa));
            float prophom2 = (1.0 * cohort_counts[2].aa/(cohort_counts[2].AA + cohort_counts[2].Aa + cohort_counts[2].aa));
            float prophom3 = (1.0 * cohort_counts[3].aa/(cohort_counts[3].AA + cohort_counts[3].Aa + cohort_counts[3].aa));
            float prophom4 = (1.0 * cohort_counts[4].aa/(cohort_counts[4].AA + cohort_counts[4].Aa + cohort_counts[4].aa));
            float prophom5 = (1.0 * cohort_counts[5].aa/(cohort_counts[5].AA + cohort_counts[5].Aa + cohort_counts[5].aa));
            fprintf(outfileprophom_handle, "%d\t%f\t%f\t%f\t%f\t%f\n", chosenFounder, prophom1, prophom2, prophom3, prophom4, prophom5); 

            // float nhom1 = cohort_counts[1].aa;
            // float nhom2 = cohort_counts[2].aa;
            // float nhom3 = cohort_counts[3].aa;
            // float nhom4 = cohort_counts[4].aa;
            // float nhom5 = cohort_counts[5].aa;
            // fprintf(outfileprophom_handle, "%d\t%f\t%f\t%f\t%f\t%f\n", chosenFounder, nhom1, nhom2, nhom3, nhom4, nhom5); 
			// reinitialize cohort_count storage
			for (i = 0; i < num_cohorts; i++) {
				init_genocount( &cohort_counts[i] );
			}
		}
	}
	
	
	fclose(outfile_handle);	
	exit(0);
}

