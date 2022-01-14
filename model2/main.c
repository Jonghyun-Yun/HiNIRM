#include "header.h"
#include "FILE/nrutil.h"
#include "FILE/stat.c"
#include "cost.c"

int main(int argc, const char * argv[])
{
	// Declare Variables
	FILE *inp, *JIN, *HUR, *OUT, *PRT, *JYW, *ASA, *IMS, *JJW;
	char buf[255], frname[255];
	int stime; long ltime;
	int ind, ite, a, b, i, j, k, l, v, *t, *tvec, accept, gcount, mcount, mutmp, *count, show1, show2;
	double num, den, un, ratio;
	double old_like_beta, new_like_beta, old_like_theta, new_like_theta;
	double update_like_samp, update_like_item, tmp_oldmu, tmp_newmu;
	double *post_ra, *post_rb, *post_fa, *post_fb, *school_a, *school_b;
	double *old_samp_distance, *new_samp_distance, *sample_samp_like;
	double *old_item_distance, *new_item_distance, *sample_item_like;
	double **sum_mu, **mu_dist, **sum_mu_dist;
	double ***sample_tau, **sum_tau, **var_tau;
	double **sample_sigma, *sum_sigma, *var_sigma;
	double ***sample_delta, **sum_delta, **var_delta;
	double ***sample_gamma, **sum_gamma, **var_gamma;
	double ***sample_varphi, **sum_varphi, **var_varphi;
	double *var_fix, *avg_fix, *var_ran, *avg_ran, *avg_beta, *var_beta;

	MM = atoi(argv[1]);

	// Set Random Seed
	ltime = time(NULL);
	stime = (unsigned int)ltime/2;
	srand(stime);
	printf("nseed = %d\n", stime);

/*	// Input Number of Thread
	# pragma omp parallel
	{
		#if defined (_OPENMP)
			k = omp_get_num_threads();
			printf("k = %d\n", k);
			srand(((unsigned int)time(NULL))^k);
		#endif
	}
*/

	// Input Parameters
	inp = fopen("DATA/parameter.txt", "r");
	if(inp == NULL) {printf("Can't open data file\n"); return 0;}
	fscanf(inp, "%d", &niter);
	fscanf(inp, "%d", &nburn);
	fscanf(inp, "%d", &thin);
	fscanf(inp, "%d", &print);
	fscanf(inp, "%d", &repeat);
	fscanf(inp, "%lf", &jump_beta);
	fscanf(inp, "%lf", &jump_theta);
	fscanf(inp, "%lf", &jump_mu);
	fscanf(inp, "%lf", &jump_W);
	fclose(inp);

	// The Number of Respondents by Schools
	ncount = ivector(1, nSCHOOL);
	inp = fopen("DATA/count.txt", "r");
	for(i = 1; i <= nSCHOOL; i++) fscanf(inp, "%d", &ncount[i]);
	fclose(inp);

	// Renovation School Systems
	t = ivector(1, nSCHOOL);
	tvec = ivector(0, nCAT);
	for(i = 0; i <= nCAT; i++) tvec[i] = 0;
	inp = fopen("DATA/renov.txt", "r");
	for(i = 1; i <= nSCHOOL; i++){
		fscanf(inp, "%d", &t[i]);
		for(j = 0; j <= nCAT; j++)
			if(t[i] == j) tvec[j]++;
	}
	fclose(inp);

	jump_Z = dvector(1, 10);
	inp = fopen("DATA/jumprule.txt", "r");
	for(i = 1; i <= 10; i++) fscanf(inp, "%lf", &jump_Z[i]);
	fclose(inp);

	jump_index = imatrix(1, nSCHOOL, 1, nITEM);
	inp = fopen("DATA/jumpitem.txt", "r");
	for(i = 1; i <= nSCHOOL; i++)
		for(j = 1; j <= nITEM; j++) fscanf(inp, "%d", &jump_index[i][j]);
	fclose(inp);

	// Declare typedef structure and set array of variables in typedef structure
	totalsize  = sizeof(SCHOOL) + sizeof(int) * (nMAX+1)*(nITEM+1);
	totalsize += sizeof(int) * (nMAX+1) + sizeof(int) * (nITEM+1);
	totalsize += sizeof(int) * (nITEM+1)*(nMAX+1)*(nMAX+1);
	totalsize += sizeof(int) * (nMAX+1)*(nITEM+1)*(nITEM+1);
	totalsize += sizeof(double) * ((nITEM+1)*2 + (nMAX+1)*2) + sizeof(double) * ((nITEM+1)*(nITEM+1)*2);
	totalsize += sizeof(double) * ((nMAX+1)*(nDIM+1)*4 + (nITEM+1)*(nDIM+1)*2);
	totalsize += sizeof(double) * (((niter-nburn)/thin+1)*(nITEM+1) + (nITEM+1)*3);
	totalsize += sizeof(double) * (((niter-nburn)/thin+1)*(nMAX+1) + (nMAX+1)*3);
	totalsize += sizeof(double) * (((niter-nburn)/thin+1)*((nMAX+1)*(nDIM+1) + (nITEM+1)*(nDIM+1)));
	totalsize += sizeof(double) * (((niter-nburn)/thin+1) + (nDIM+1));
	totalsize += sizeof(double) * ((nMAX+1)*(nDIM+1)*2 + (nITEM+1)*(nDIM+1)*2 + (nMAX+1) + (nITEM+1));
	totalsize += sizeof(double) * ((nITEM+1)*(nITEM+1)*3);
	SCHOOL = (YEWON *)malloc(totalsize * (nSCHOOL+1));
	for(k = 0; k <= nSCHOOL; k++){
		SCHOOL[k].cbsize = totalsize;
		SCHOOL[k].dataset = (int**)malloc(sizeof(int*)*(nMAX+1));
		SCHOOL[k].count_samp = (int*)malloc(sizeof(int*)*(nMAX+1));
		SCHOOL[k].count_item = (int*)malloc(sizeof(int*)*(nITEM+1));
		SCHOOL[k].Y = (int***)malloc(sizeof(int**)*(nITEM+1));
		SCHOOL[k].U = (int***)malloc(sizeof(int**)*(nMAX+1));
		SCHOOL[k].oldbeta = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].newbeta = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].oldtheta = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].newtheta = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].old_Zsamp = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].new_Zsamp = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].old_Zmean = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].new_Zmean = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].old_Zitem = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].new_Zitem = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].mean_Z = (double*)malloc(sizeof(double)*(nDIM+1));
		SCHOOL[k].sample_beta  = (double**)malloc(sizeof(double*)*((niter-nburn)/thin+1));
		SCHOOL[k].sample_theta = (double**)malloc(sizeof(double*)*((niter-nburn)/thin+1));
		SCHOOL[k].sample_sigma = (double*)malloc(sizeof(double)*((niter-nburn)/thin+1));
		SCHOOL[k].sum_beta = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].var_beta = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].acc_beta = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].sum_theta = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].var_theta = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].acc_theta = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].sample_Zsamp = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
		SCHOOL[k].sample_Zitem = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
		SCHOOL[k].sample_item_mat = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].sum_Zsamp = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].var_Zsamp = (double**)malloc(sizeof(double*)*(nMAX+1));
		SCHOOL[k].acc_Zsamp = (double*)malloc(sizeof(double)*(nMAX+1));
		SCHOOL[k].sum_Zitem = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].var_Zitem = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].acc_Zitem = (double*)malloc(sizeof(double)*(nITEM+1));
		SCHOOL[k].old_item_mat = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].new_item_mat = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].sum_item_mat = (double**)malloc(sizeof(double*)*(nITEM+1));
		SCHOOL[k].var_item_mat = (double**)malloc(sizeof(double*)*(nITEM+1));
		for(i = 0; i <= nMAX; i++) SCHOOL[k].dataset[i] = (int*)malloc(sizeof(int)*(nITEM+1));
		for(i = 0; i <= nITEM; i++){
			SCHOOL[k].Y[i] = (int**)malloc(sizeof(int*)*(nMAX+1));
			for(a = 0; a <= nMAX; a++) SCHOOL[k].Y[i][a] = (int*)malloc(sizeof(int)*(nMAX+1));
		}
		for(i = 0; i <= nMAX; i++){
			SCHOOL[k].U[i] = (int**)malloc(sizeof(int*)*(nITEM+1));
			for(a = 0; a <= nITEM; a++) SCHOOL[k].U[i][a] = (int*)malloc(sizeof(int)*(nITEM+1));
		}
		for(i = 0; i <= nMAX; i++){
			SCHOOL[k].old_Zsamp[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].new_Zsamp[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].old_Zmean[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].new_Zmean[i] = (double*)malloc(sizeof(double)*(nDIM+1));
		}
		for(i = 0; i <= nITEM; i++){
			SCHOOL[k].old_Zitem[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].new_Zitem[i] = (double*)malloc(sizeof(double)*(nDIM+1));
		}
		for(i = 0; i <= (niter-nburn)/thin; i++){
			SCHOOL[k].sample_beta[i]  = (double*)malloc(sizeof(double)*(nITEM+1));
			SCHOOL[k].sample_theta[i] = (double*)malloc(sizeof(double)*(nMAX+1));
			SCHOOL[k].sample_Zsamp[i] = (double**)malloc(sizeof(double*)*(nMAX+1));
			SCHOOL[k].sample_Zitem[i] = (double**)malloc(sizeof(double*)*(nITEM+1));
			for(j = 0; j <= nMAX; j++)  SCHOOL[k].sample_Zsamp[i][j] = (double*)malloc(sizeof(double)*(nDIM+1));
			for(j = 0; j <= nITEM; j++) SCHOOL[k].sample_Zitem[i][j] = (double*)malloc(sizeof(double)*(nDIM+1));
		}
		for(i = 0; i <= nMAX; i++){
			SCHOOL[k].sum_Zsamp[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].var_Zsamp[i] = (double*)malloc(sizeof(double)*(nDIM+1));
		}
		for(i = 0; i <= nITEM; i++){
			SCHOOL[k].sum_Zitem[i] = (double*)malloc(sizeof(double)*(nDIM+1));
			SCHOOL[k].var_Zitem[i] = (double*)malloc(sizeof(double)*(nDIM+1));
		}
		for(i = 0; i <= nITEM; i++){
			SCHOOL[k].sample_item_mat[i] = (double*)malloc(sizeof(double)*(nITEM+1));
			SCHOOL[k].old_item_mat[i] = (double*)malloc(sizeof(double)*(nITEM+1));
			SCHOOL[k].new_item_mat[i] = (double*)malloc(sizeof(double)*(nITEM+1));
			SCHOOL[k].sum_item_mat[i] = (double*)malloc(sizeof(double)*(nITEM+1));
			SCHOOL[k].var_item_mat[i] = (double*)malloc(sizeof(double)*(nITEM+1));
		}
		printf("MEMORY SETTING: %.2d\n", k);
	}

	count = ivector(1, nSCHOOL);

	oldmu = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nSCHOOL);
	olddelta = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	oldsigma = dvector(1, nSCHOOL);
	oldtau = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	oldgamma = dmatrix(1, nITEM, 0, nCAT);
	oldvarphi = dmatrix(1, nITEM, 0, nCAT);

	sample_sigma = dmatrix(1, (niter-nburn)/thin, 1, nSCHOOL);
	sample_delta = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
	sample_tau = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
	sample_gamma = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
	sample_varphi = (double***)malloc(sizeof(double**)*((niter-nburn)/thin+1));
	for(i = 0; i <= (niter-nburn)/thin; i++){
		sample_delta[i] = (double**)malloc(sizeof(double*)*(nITEM*(nITEM-1)/2+1));
		sample_tau[i] = (double**)malloc(sizeof(double*)*(nITEM*(nITEM-1)/2+1));
		sample_gamma[i] = (double**)malloc(sizeof(double*)*(nITEM+1));
		sample_varphi[i] = (double**)malloc(sizeof(double*)*(nITEM+1));
		for(j = 0; j <= nITEM * (nITEM - 1) / 2; j++){
			sample_delta[i][j] = (double*)malloc(sizeof(double)*(nCAT+1));
			sample_tau[i][j] = (double*)malloc(sizeof(double)*(nCAT+1));
		}
		for(j = 0; j <= nITEM; j++){
			sample_gamma[i][j] = (double*)malloc(sizeof(double)*(nCAT+1));
			sample_varphi[i][j] = (double*)malloc(sizeof(double)*(nCAT+1));
		}
	}
	sum_mu = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nSCHOOL);
	sum_tau = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	var_tau = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	sum_sigma = dvector(1, nSCHOOL);
	var_sigma = dvector(1, nSCHOOL);
	sum_delta = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	var_delta = dmatrix(1, nITEM * (nITEM - 1) / 2, 0, nCAT);
	sum_gamma = dmatrix(1, nITEM, 0, nCAT);
	var_gamma = dmatrix(1, nITEM, 0, nCAT);
	sum_varphi = dmatrix(1, nITEM, 0, nCAT);
	var_varphi = dmatrix(1, nITEM, 0, nCAT);

	post_ra = dvector(1, nSCHOOL);
	post_rb = dvector(1, nSCHOOL);
	post_fa = dvector(0, nCAT);
	post_fb = dvector(0, nCAT);
	school_a = dvector(0, nCAT);
	school_b = dvector(0, nCAT);

	var_fix = dvector(0, nCAT);
	avg_fix = dvector(0, nCAT);
	var_ran = dvector(1, nSCHOOL);
	avg_ran = dvector(1, nSCHOOL);
	avg_beta = dvector(0, nCAT);
	var_beta = dvector(0, nCAT);

	mu_dist = dmatrix(1, nSCHOOL, 1, nSCHOOL);
	sum_mu_dist = dmatrix(1, nSCHOOL, 1, nSCHOOL);

	frname[0]  = 'D'; frname[1]  = 'A'; frname[2]  = 'T'; frname[3]  = 'A'; frname[4]  = '/';
	frname[5]  = 'i'; frname[6]  = 't'; frname[7]  = 'e'; frname[8]  = 'm';
	frname[11] = '.'; frname[12] = 't'; frname[13] = 'x'; frname[14] = 't'; frname[15] = '\0';

	for(k = 0; k <= nSCHOOL; k++){
		for(i = 0; i <= nMAX; i++) SCHOOL[k].count_samp[i] = 0;
		for(i = 0; i <= nITEM; i++) SCHOOL[k].count_item[i] = 0;
		for(i = 0; i <= nMAX; i++)
			for(j = 0; j <= nITEM; j++) SCHOOL[k].dataset[i][j] = 0;
		for(i = 0; i <= nITEM; i++) SCHOOL[k].oldbeta[i] = SCHOOL[k].newbeta[i] = 0.0;
		for(i = 0; i <= nMAX; i++) SCHOOL[k].oldtheta[i] = SCHOOL[k].newtheta[i] = 0.0;
		for(i = 0; i <= nITEM; i++)
			for(j = 0; j <= nITEM; j++) SCHOOL[k].old_item_mat[i][j] = SCHOOL[k].new_item_mat[i][j] = 0.0;
		for(i = 0; i <= nITEM; i++)
			for(a = 0; a <= nMAX; a++)
				for(b = 0; b <= nMAX; b++) SCHOOL[k].Y[i][a][b] = 0;
		for(i = 0; i <= nMAX; i++)
			for(a = 0; a <= nITEM; a++)
				for(b = 0; b <= nITEM; b++) SCHOOL[k].U[i][a][b] = 0;
		for(i = 0; i <= nMAX; i++)
			for(j = 0; j <= nDIM; j++) SCHOOL[k].old_Zsamp[i][j] = SCHOOL[k].new_Zsamp[i][j] = SCHOOL[k].old_Zmean[i][j] = SCHOOL[k].new_Zmean[i][j] = 0.0;
		for(i = 0; i <= nITEM; i++)
			for(j = 0; j <= nDIM; j++) SCHOOL[k].old_Zitem[i][j] = SCHOOL[k].new_Zitem[i][j] = 0.0;
		for(i = 0; i <= (niter-nburn)/thin; i++){
			SCHOOL[k].sample_sigma[i] = 0.0;
			for(j = 0; j <= nITEM; j++) SCHOOL[k].sample_beta[i][j] = 0.0;
			for(j = 0; j <= nMAX; j++) SCHOOL[k].sample_theta[i][j] = 0.0;
			for(a = 0; a <= nMAX; a++)
				for(b = 0; b <= nDIM; b++) SCHOOL[k].sample_Zsamp[i][a][b] = 0.0;
			for(a = 0; a <= nITEM; a++)
				for(b = 0; b <= nDIM; b++) SCHOOL[k].sample_Zitem[i][a][b] = 0.0;
		}
		SCHOOL[k].oldsigma = SCHOOL[k].sum_sigma = SCHOOL[k].var_sigma = 0.0;
		for(i = 0; i <= nDIM; i++) SCHOOL[k].mean_Z[i] = 0.0;
		for(i = 0; i <= nITEM; i++) SCHOOL[k].var_beta[i] = SCHOOL[k].sum_beta[i] = SCHOOL[k].acc_beta[i] = 0.0;
		for(i = 0; i <= nMAX; i++) SCHOOL[k].var_theta[i] = SCHOOL[k].sum_theta[i] = SCHOOL[k].acc_theta[i] = 0.0;
		for(i = 0; i <= nMAX; i++)
			for(j = 0; j <= nDIM; j++) SCHOOL[k].sum_Zsamp[i][j] = SCHOOL[k].var_Zsamp[i][j] = 0.0;
		for(i = 0; i <= nITEM; i++)
			for(j = 0; j <= nDIM; j++) SCHOOL[k].sum_Zitem[i][j] = SCHOOL[k].var_Zitem[i][j] = 0.0;
		for(i = 0; i <= nITEM; i++)
			for(j = 0; j <= nITEM; j++) SCHOOL[k].sample_item_mat[i][j] =  SCHOOL[k].sum_item_mat[i][j] = SCHOOL[k].var_item_mat[i][j] = 0.0;
		for(i = 0; i <= nITEM; i++) SCHOOL[k].acc_Zitem[i] = 0.0;
		for(i = 0; i <= nMAX; i++) SCHOOL[k].acc_Zsamp[i] = 0.0;
		if(k != 0) count[k] = 0;

		if(k != 0){
			if(k < 10){frname[9] = (char)(48); frname[10] = (char)(k + 48);}
			else{frname[9] = (char)(k/10 + 48); frname[10] = (char)(k%10 + 48);}

			inp = fopen(frname, "r"); printf("Currently Reading %s\n", frname);
			if(inp == NULL) {printf("Cannot open data file\n"); return 0;}
			for(i = 1; i <= ncount[k]; i++)
				for(j = 1; j <= nITEM; j++){
					fscanf(inp, "%d", &SCHOOL[k].dataset[i][j]);
					SCHOOL[k].count_samp[i] += SCHOOL[k].dataset[i][j];
					SCHOOL[k].count_item[j] += SCHOOL[k].dataset[i][j];
				}
			fclose(inp);

			printf("%.2d\n", k);
			for(i = 1; i <= ncount[k]; i++){
				for(j = 1; j <= nITEM; j++) printf("%d ", SCHOOL[k].dataset[i][j]);
				printf("\n");
			}

			for(i = 1; i <= nITEM; i++)
				for(a = 2; a <= ncount[k]; a++)
					for(b = 1; b < a; b++){
						SCHOOL[k].Y[i][a][b] = SCHOOL[k].dataset[a][i] * SCHOOL[k].dataset[b][i];
						SCHOOL[k].Y[i][b][a] = SCHOOL[k].Y[i][a][b];
					}

			for(a = 1; a <= ncount[k]; a++)
				for(i = 2; i <= nITEM; i++)
					for(j = 1; j < i; j++){
						SCHOOL[k].U[a][i][j] = SCHOOL[k].dataset[a][i] * SCHOOL[k].dataset[a][j];
						SCHOOL[k].U[a][j][i] = SCHOOL[k].U[a][i][j];
					}
		}
		printf("INITIALIZATION AND DATA LOADING: %.2d\n", k);
	}

	// Declare Additional Variables
	sample_samp_like = dvector(1, nMAX);
	old_samp_distance = dvector(1, nMAX);
	new_samp_distance = dvector(1, nMAX);

	sample_item_like = dvector(1, nITEM);
	old_item_distance = dvector(1, nITEM);
	new_item_distance = dvector(1, nITEM);
	pr_var_Z = sqrt(2.0);

	for(v = 0; v < repeat; v++){
		// Initialize Variables
		for(k = 1; k <= nSCHOOL; k++){
			for(i = 1; i <= nITEM; i++) SCHOOL[k].oldbeta[i] = SCHOOL[k].newbeta[i] = 0.0;
			for(i = 1; i <= ncount[k]; i++) SCHOOL[k].oldtheta[i] = SCHOOL[k].newtheta[i] = 0.0;
			for(i = 1; i <= ncount[k]; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].old_Zsamp[i][j] = SCHOOL[k].new_Zsamp[i][j] = SCHOOL[k].old_Zmean[i][j] = SCHOOL[k].new_Zmean[i][j] = 0.0;;
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].old_Zitem[i][j] = SCHOOL[k].new_Zitem[i][j] = 0.0;
			for(i = 1; i <= (niter-nburn)/thin; i++){
				SCHOOL[k].sample_sigma[i] = 0.0;
				for(j = 1; j <= nITEM; j++) SCHOOL[k].sample_beta[i][j] = 0.0;
				for(j = 1; j <= ncount[k]; j++) SCHOOL[k].sample_theta[i][j] = 0.0;
				for(a = 1; a <= ncount[k]; a++)
					for(b = 1; b <= nDIM; b++) SCHOOL[k].sample_Zsamp[i][a][b] = 0.0;
				for(a = 1; a <= nITEM; a++)
					for(b = 1; b <= nDIM; b++) SCHOOL[k].sample_Zitem[i][a][b] = 0.0;
			}
			for(i = 1; i <= nITEM; i++) SCHOOL[k].var_beta[i] = SCHOOL[k].sum_beta[i] = SCHOOL[k].acc_beta[i] = 0.0;
			for(i = 1; i <= ncount[k]; i++) SCHOOL[k].var_theta[i] = SCHOOL[k].sum_theta[i] = SCHOOL[k].acc_theta[i] = 0.0;
			for(i = 1; i <= ncount[k]; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].sum_Zsamp[i][j] = SCHOOL[k].var_Zsamp[i][j] = 0.0;
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].sum_Zitem[i][j] = SCHOOL[k].var_Zitem[i][j] = 0.0;
			for(i = 1; i <= ncount[k]; i++) SCHOOL[k].acc_Zsamp[i] = 0.0;
			for(i = 0; i <= nITEM; i++) SCHOOL[k].acc_Zitem[i] = 0.0;
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nITEM; j++){
					SCHOOL[k].sample_item_mat[i][j] = 0.0;
					SCHOOL[k].old_item_mat[i][j] = SCHOOL[k].new_item_mat[i][j] = 0.0;
					SCHOOL[k].sum_item_mat[i][j] = SCHOOL[k].var_item_mat[i][j] = 0.0;
				}
			for(i = 0; i <= nDIM; i++) SCHOOL[k].mean_Z[i] = 0.0;
			SCHOOL[k].oldsigma = SCHOOL[k].sum_sigma = SCHOOL[k].var_sigma = 0.0;
			count[k] = 0;
		}

		for(k = 1; k <= nSCHOOL; k++){
			oldsigma[k] = 0.0;
			sum_sigma[k] = var_sigma[k] = 0.0;
			for(j = 1; j <= (niter-nburn)/thin; j++) sample_sigma[j][k] = 0.0;
		}
		for(k = 0; k <= nCAT; k++){
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++){
				olddelta[i][k] = oldtau[i][k] = 0.0;
				sum_delta[i][k] = var_delta[i][k] = 0.0;
				sum_tau[i][k] = var_tau[i][k] = 0.0;
				for(j = 1; j <= (niter-nburn)/thin; j++) sample_tau[j][i][k] = sample_delta[j][i][k] = 0.0;
			}
			for(i = 1; i <= nITEM; i++){
				oldgamma[i][k] = oldvarphi[i][k] = 0.0;
				sum_gamma[i][k] = var_gamma[i][k] = 0.0;
				sum_varphi[i][k] = var_varphi[i][k] = 0.0;
				for(j = 1; j <= (niter-nburn)/thin; j++) sample_gamma[j][i][k] = sample_varphi[j][i][k] = 0.0;
			}
		}
		for(k = 1; k <= nSCHOOL; k++)
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) oldmu[i][k] = sum_mu[i][k] = 0.0;
		for(i = 1; i <= nSCHOOL; i++)
			for(j = 1; j <= nSCHOOL; j++) sum_mu_dist[i][j] = 0.0;

		// Generate Initial Values for beta, Z, sigma
		for(i = 1; i <= nSCHOOL; i++) oldsigma[i] = 100.0;
		for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++){
			for(j = 0; j <= nCAT; j++){
				olddelta[i][j] = -1.5 + 3.0 * rand() / RAND_MAX;
				oldtau[i][j] = 100.0;
			}
			for(j = 1; j <= nSCHOOL; j++) oldmu[i][j] = -1.5 + 3.0 * rand() / RAND_MAX;
		}
		for(i = 1; i <= nITEM; i++)
			for(j = 0; j <= nCAT; j++){
				oldgamma[i][j] = -1.5 + 3.0 * rand() / RAND_MAX;
				oldvarphi[i][j] = 100.0;
			}

		for(k = 1; k <= nSCHOOL; k++){
			SCHOOL[k].oldsigma = 0.05 * 0.05;
			for(i = 1; i <= nITEM; i++) SCHOOL[k].oldbeta[i] = -1.5 + 3.0 * rand() / RAND_MAX;
			for(i = 1; i <= ncount[k]; i++) SCHOOL[k].oldtheta[i] = -1.5 + 3.0 * rand() / RAND_MAX;
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].old_Zitem[i][j] = SCHOOL[k].new_Zitem[i][j] = -1.5 + 3.0 * rand() / RAND_MAX;
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++)
					for(a = 1; a <= ncount[k]; a++)
						if(SCHOOL[k].dataset[a][i] == 1) SCHOOL[k].old_Zmean[a][j] += SCHOOL[k].old_Zitem[i][j] / (SCHOOL[k].count_samp[a] * 1.0);
			for(i = 1; i <= ncount[k]; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].new_Zmean[i][j] = SCHOOL[k].old_Zmean[i][j];
			for(i = 1; i <= ncount[k]; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[k].new_Zsamp[i][j] = SCHOOL[k].old_Zsamp[i][j] = SCHOOL[k].old_Zmean[i][j] + sqrt(SCHOOL[k].oldsigma) * gasdev();
			for(i = 2; i <= nITEM; i++)
				for(j = 1; j < i; j++){
					for(l = 1; l <= nDIM; l++) SCHOOL[k].old_item_mat[i][j] += pow((SCHOOL[k].old_Zitem[i][l] - SCHOOL[k].old_Zitem[j][l]), 2.0);
					SCHOOL[k].old_item_mat[i][j] = sqrt(SCHOOL[k].old_item_mat[i][j]);
					SCHOOL[k].old_item_mat[j][i] = SCHOOL[k].old_item_mat[i][j];
				}
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nITEM; j++) SCHOOL[k].new_item_mat[i][j] = SCHOOL[k].old_item_mat[i][j];
		}
		pr_var_Z = log(sqrt(2.0));

		// MCMC Implementation for Parameter Estimation
		frname[0] = 'R'; frname[1] = 'E'; frname[2] = 'S'; frname[3] = 'U'; frname[4] = 'L'; frname[5] = 'T'; frname[6] = '/';
		frname[7] = 's'; frname[8] = 'i'; frname[9] = 'm'; frname[10] = '_'; frname[12] = (char)(48+MM); frname[13] = '.'; frname[14] = 'l'; frname[15] = 'o'; frname[16] = 'g'; frname[17] = '\0';
		frname[11] = 's'; HUR = fopen(frname, "a");
		frname[11] = 'l'; JYW = fopen(frname, "a");
		frname[11] = 'u'; OUT = fopen(frname, "a");
		frname[11] = 'g'; JIN = fopen(frname, "a");
		frname[11] = 'p'; PRT = fopen(frname, "a");
		frname[11] = 'a'; JJW = fopen(frname, "a");
		gcount = mcount = 0;
		for(iter = 1; iter <= niter; iter++){
			for(a = 1; a <= nSCHOOL; a++){

				for(i = 1; i <= nITEM; i++){
					//#pragma omp parallel for private(j, k) default(shared)
					for(j = 1; j <= nDIM; j++){
						SCHOOL[a].new_Zitem[i][j] = SCHOOL[a].old_Zitem[i][j] + jump_Z[jump_index[a][i]] * gasdev();
						for(k = 1; k <= ncount[a]; k++)
							if(SCHOOL[a].dataset[k][i] == 1){
								SCHOOL[a].new_Zmean[k][j] -= SCHOOL[a].old_Zitem[i][j] / (SCHOOL[a].count_samp[k] * 1.0);
								SCHOOL[a].new_Zmean[k][j] += SCHOOL[a].new_Zitem[i][j] / (SCHOOL[a].count_samp[k] * 1.0);
							}
					}

					for(ind = 1; ind <= nITEM; ind++) sample_item_like[ind] = old_item_distance[ind] = new_item_distance[ind] = 0.0;
					//#pragma omp parallel for private(ind, k, l) default(shared)
					for(ind = 1; ind <= nITEM; ind++)
						if(ind != i){
							for(l = 1; l <= nDIM; l++){
								old_item_distance[ind] += pow((SCHOOL[a].old_Zitem[ind][l] - SCHOOL[a].old_Zitem[i][l]), 2.0);
								new_item_distance[ind] += pow((SCHOOL[a].new_Zitem[ind][l] - SCHOOL[a].new_Zitem[i][l]), 2.0);
							}
							old_item_distance[ind] = sqrt(old_item_distance[ind]);
							new_item_distance[ind] = sqrt(new_item_distance[ind]);
							SCHOOL[a].new_item_mat[ind][i] = new_item_distance[ind];
							SCHOOL[a].new_item_mat[i][ind] = SCHOOL[a].new_item_mat[ind][i];
							SCHOOL[a].old_item_mat[ind][i] = old_item_distance[ind];
							SCHOOL[a].old_item_mat[i][ind] = SCHOOL[a].old_item_mat[ind][i];
							for(k = 1; k <= ncount[a]; k++){
								if(SCHOOL[a].U[k][ind][i] == 1){
									sample_item_like[ind] -= -log(1.0 + exp(-(SCHOOL[a].oldtheta[k] - old_item_distance[ind])));
									sample_item_like[ind] += -log(1.0 + exp(-(SCHOOL[a].oldtheta[k] - new_item_distance[ind])));
								}
								else{
									sample_item_like[ind] -= -log(1.0 + exp(SCHOOL[a].oldtheta[k] - old_item_distance[ind]));
									sample_item_like[ind] += -log(1.0 + exp(SCHOOL[a].oldtheta[k] - new_item_distance[ind]));
								}
							}
						}
					update_like_item = 0.0;
					for(ind = 1; ind <= nITEM; ind++) update_like_item += sample_item_like[ind];

					num = den = 0.0;
					for(j = 2; j <= nITEM; j++)
						for(k = 1; k < j; k++){
							if(SCHOOL[a].new_item_mat[j][k] > 0.0001) num += dlognorm(log(SCHOOL[a].new_item_mat[j][k]), olddelta[((j-1)*(j-2)/2+k)][t[a]], sqrt(oldtau[((j-1)*(j-2)/2+k)][t[a]]));
							else num += dlognorm(log(0.0001), olddelta[((j-1)*(j-2)/2+k)][t[a]], sqrt(oldtau[((j-1)*(j-2)/2+k)][t[a]]));
							if(SCHOOL[a].old_item_mat[j][k] > 0.0001) den += dlognorm(log(SCHOOL[a].old_item_mat[j][k]), olddelta[((j-1)*(j-2)/2+k)][t[a]], sqrt(oldtau[((j-1)*(j-2)/2+k)][t[a]]));
							else den += dlognorm(log(0.0001), olddelta[((j-1)*(j-2)/2+k)][t[a]], sqrt(oldtau[((j-1)*(j-2)/2+k)][t[a]]));
						}

					ratio = update_like_item + (num - den);
					//printf("SCHOOL-%.2d, ITEM-%.2d: Num-%.3f, Den-%.3f\n", a, i, num, den);

					if(ratio > 0.0) accept = 1;
					else{
						un = rand() * 1.0 / RAND_MAX;
						if(log(un) < ratio) accept = 1;
						else accept = 0;
					}

					if(accept == 1){
						for(j = 1; j <= nDIM; j++){
							SCHOOL[a].old_Zitem[i][j] = SCHOOL[a].new_Zitem[i][j];
							for(k = 1; k <= ncount[a]; k++)
								if(SCHOOL[a].dataset[k][i] == 1) SCHOOL[a].old_Zmean[k][j] = SCHOOL[a].new_Zmean[k][j];
						}
						SCHOOL[a].acc_Zitem[i] += 1.0 / niter;
						for(j = 1; j <= nITEM; j++)
							for(k = 1; k <= nITEM; k++) SCHOOL[a].old_item_mat[j][k] = SCHOOL[a].new_item_mat[j][k];
					}
					else{
						for(j = 1; j <= nDIM; j++){
							SCHOOL[a].new_Zitem[i][j] = SCHOOL[a].old_Zitem[i][j];
							for(k = 1; k <= ncount[a]; k++)
								if(SCHOOL[a].dataset[k][i] == 1) SCHOOL[a].new_Zmean[k][j] = SCHOOL[a].old_Zmean[k][j];
						}
						for(j = 1; j <= nITEM; j++)
							for(k = 1; k <= nITEM; k++) SCHOOL[a].new_item_mat[j][k] = SCHOOL[a].old_item_mat[j][k];
					}
				}

				for(i = 1; i <= ncount[a]; i++){
					for(j = 1; j <= nDIM; j++) SCHOOL[a].new_Zsamp[i][j] = SCHOOL[a].old_Zsamp[i][j] + jump_W * gasdev();

					for(ind = 1; ind <= ncount[a]; ind++) sample_samp_like[ind] = old_samp_distance[ind] = new_samp_distance[ind] = 0.0;
					//#pragma omp parallel for private(ind, k, l) default(shared)
					for(ind = 1; ind <= ncount[a]; ind++)
						if(ind != i){
							for(l = 1; l <= nDIM; l++){
								old_samp_distance[ind] += pow((SCHOOL[a].old_Zsamp[ind][l] - SCHOOL[a].old_Zsamp[i][l]), 2.0);
								new_samp_distance[ind] += pow((SCHOOL[a].old_Zsamp[ind][l] - SCHOOL[a].new_Zsamp[i][l]), 2.0);
							}
							old_samp_distance[ind] = sqrt(old_samp_distance[ind]);
							new_samp_distance[ind] = sqrt(new_samp_distance[ind]);
							for(k = 1; k <= nITEM; k++){
								if(SCHOOL[a].Y[k][ind][i] == 1){
									sample_samp_like[ind] -= -log(1.0 + exp(-(SCHOOL[a].oldbeta[k] - old_samp_distance[ind])));
									sample_samp_like[ind] += -log(1.0 + exp(-(SCHOOL[a].oldbeta[k] - new_samp_distance[ind])));
								}
								else{
									sample_samp_like[ind] -= -log(1.0 + exp(SCHOOL[a].oldbeta[k] - old_samp_distance[ind]));
									sample_samp_like[ind] += -log(1.0 + exp(SCHOOL[a].oldbeta[k] - new_samp_distance[ind]));
								}
							}
						}
					update_like_samp = 0.0;
					for(ind = 1; ind <= ncount[a]; ind++) update_like_samp += sample_samp_like[ind];
					//printf("SCHOOL-%.2d, PERSON-%.2d: LIKELIHOOD_PERSON-%.3f\n", a, i, update_like_samp);

					num = den = 0.0;
					//printf("SCHOOL-%.2d, PERSON-%.2d: Num-%.3f, Den-%.3f\n", a, i, num, den);
					for(j = 1; j <= nDIM; j++){
						num += dlognorm(SCHOOL[a].new_Zsamp[i][j], SCHOOL[a].old_Zmean[i][j], sqrt(SCHOOL[a].oldsigma));
						den += dlognorm(SCHOOL[a].old_Zsamp[i][j], SCHOOL[a].old_Zmean[i][j], sqrt(SCHOOL[a].oldsigma));
					}
					ratio = update_like_samp + (num - den);
					//printf("SCHOOL-%.2d, PERSON-%.2d: Num-%.3f, Den-%.3f\n", a, i, num, den);

					if(ratio > 0.0) accept = 1;
					else{
						un = rand() * 1.0 / RAND_MAX;
						if(log(un) < ratio) accept = 1;
						else accept = 0;
					}

					if(accept == 1){
						for(j = 1; j <= nDIM; j++) SCHOOL[a].old_Zsamp[i][j] = SCHOOL[a].new_Zsamp[i][j];
						SCHOOL[a].acc_Zsamp[i] += 1.0 / niter;
					}
					else{
						for(j = 1; j <= nDIM; j++) SCHOOL[a].new_Zsamp[i][j] = SCHOOL[a].old_Zsamp[i][j];
					}
				}

				SCHOOL[a].post_a = prior_a; SCHOOL[a].post_b = prior_b;
				for(i = 1; i <= ncount[a]; i++)
					for(j = 1; j <= nDIM; j++){
						SCHOOL[a].post_a += 0.5;
						SCHOOL[a].post_b += 0.5 * (SCHOOL[a].old_Zsamp[i][j] - SCHOOL[a].old_Zmean[i][j]) * (SCHOOL[a].old_Zsamp[i][j] - SCHOOL[a].old_Zmean[i][j]);
					}
				SCHOOL[a].oldsigma = 1.0 / Rgamma(SCHOOL[a].post_a, SCHOOL[a].post_b);

				// 2. Update $\beta_i$ from the proposal distribution $\phi_2(\cdot)$
				//#pragma omp parallel for private(i, j, k, old_like_beta, new_like_beta, num, den, accept, ratio, un) default(shared)
				for(i = 1; i <= nITEM; i++){
					old_like_beta = cost_beta(i, SCHOOL[a].oldbeta[i], a);
					SCHOOL[a].newbeta[i] = SCHOOL[a].oldbeta[i] + jump_beta * gasdev();
					if(fabs(SCHOOL[a].newbeta[i]) < 7.0){
						new_like_beta = cost_beta(i, SCHOOL[a].newbeta[i], a);
						num = new_like_beta; den = old_like_beta;

						num += dlognorm(SCHOOL[a].oldbeta[i], oldgamma[i][t[a]], sqrt(oldvarphi[i][t[a]]));
						den += dlognorm(SCHOOL[a].newbeta[i], oldgamma[i][t[a]], sqrt(oldvarphi[i][t[a]]));
						ratio = num - den;

						if(ratio > 0.0) accept = 1;
						else{
							un = rand() * 1.0 / RAND_MAX;
							if(log(un) < ratio) accept = 1;
							else accept = 0;
						}
					}
					else accept = 0;

					if(accept == 1){
						SCHOOL[a].oldbeta[i] = SCHOOL[a].newbeta[i];
						SCHOOL[a].acc_beta[i] += 1.0 / niter;
					}
					else SCHOOL[a].newbeta[i] = SCHOOL[a].oldbeta[i];
				}

				//#pragma omp parallel for private(i, old_like_theta, new_like_theta, num, den, accept, ratio, un) default(shared)
				for(i = 1; i <= ncount[a]; i++){
					old_like_theta = cost_theta(i, SCHOOL[a].oldtheta[i], a);
					SCHOOL[a].newtheta[i] = SCHOOL[a].oldtheta[i] + jump_theta * gasdev();
					new_like_theta = cost_theta(i, SCHOOL[a].newtheta[i], a);

					num = dlognorm(SCHOOL[a].newtheta[i], pr_mean_theta, pr_var_theta) + new_like_theta;
					den = dlognorm(SCHOOL[a].oldtheta[i], pr_mean_theta, pr_var_theta) + old_like_theta;
					ratio = num - den;

					if(ratio > 0.0) accept = 1;
					else{
						un = rand() * 1.0 / RAND_MAX;
						if(log(un) < ratio) accept = 1;
						else accept = 0;
					}

					if(accept == 1){
						SCHOOL[a].oldtheta[i] = SCHOOL[a].newtheta[i];
						SCHOOL[a].acc_theta[i] += 1.0 / niter;
					}
					else SCHOOL[a].newtheta[i] = SCHOOL[a].oldtheta[i];
				}

				// Save MCMC Results to Files and Repository Variables
				if(iter > nburn && iter % thin == 0){
					count[a]++;
					for(i = 1; i <= ncount[a]; i++)
						for(j = 1; j <= nDIM; j++) SCHOOL[a].sample_Zsamp[count[a]][i][j] = SCHOOL[a].old_Zsamp[i][j];
					for(i = 1; i <= nITEM; i++)
						for(j = 1; j <= nDIM; j++) SCHOOL[a].sample_Zitem[count[a]][i][j] = SCHOOL[a].old_Zitem[i][j];
					for(i = 1; i <= nITEM; i++) SCHOOL[a].sample_beta[count[a]][i] = SCHOOL[a].oldbeta[i];
					for(i = 1; i <= ncount[a]; i++) SCHOOL[a].sample_theta[count[a]][i] = SCHOOL[a].oldtheta[i];
					SCHOOL[a].sample_sigma[count[a]] = SCHOOL[a].oldsigma;
				}

				// Print MCMC Results to Screen
				if(iter % print == 0){
					printf("%.5d-BETA%.2d  ", iter, a); for(i = 1; i <= nITEM; i++) printf("% .4f ", SCHOOL[a].oldbeta[i]); printf("%.4f\n", SCHOOL[a].oldsigma);
				}
			}

			for(i = 1; i <= nITEM; i++){
				for(k = 0; k <= nCAT; k++){school_a[k] = prior_a; school_b[k] = prior_b;}
				for(j = 1; j <= nSCHOOL; j++)
					for(k = 0; k <= nCAT; k++)
						if(t[j] == k){
							school_a[k] += 0.5;
							school_b[k] += 0.5 * (SCHOOL[j].oldbeta[i] - oldgamma[i][k]) * (SCHOOL[j].oldbeta[i] - oldgamma[i][k]);
						}
				for(k = 0; k <= nCAT; k++) oldvarphi[i][k] = 1.0 / Rgamma(school_a[k], school_b[k]);

				for(k = 0; k <= nCAT; k++){
					var_beta[k] = 1.0 / (1.0 / pr_var_gamma + tvec[k] / oldvarphi[i][k]);
					avg_beta[k] = 0.0;
					for(j = 1; j <= nSCHOOL; j++)
					 	if(t[j] == k) avg_beta[k] += SCHOOL[j].oldbeta[i] / tvec[k];
					avg_beta[k] *= var_beta[k] * (tvec[k] / oldvarphi[i][k]);
					oldgamma[i][k] = avg_beta[k] + sqrt(var_beta[k]) * gasdev();
				}
			}
			if(iter % print == 0){
				for(i = 1; i <= nITEM; i++){
					printf("ITERATION %.5d, GAMMA, VARPHI - ITEM%.2d: ", iter, i);
					for(j = 0; j <= nCAT; j++) printf("% .4f ", oldgamma[i][j]);
					for(j = 0; j <= nCAT; j++) printf("% .4f ", oldvarphi[i][j]);
					printf("\n");
				}
			}

			for(i = 2; i <= nITEM; i++)
				for(j = 1; j < i; j++){
					for(k = 0; k <= nCAT; k++){post_fa[k] = prior_a; post_fb[k] = prior_b;}
					//#pragma omp parallel for private(k, l) default(shared)
					for(k = 1; k <= nSCHOOL; k++)
						for(l = 0; l <= nCAT; l++)
							if(t[k] == l){
								post_fa[l] += 0.5;
								if(SCHOOL[k].old_item_mat[i][j] > 0.0001) post_fb[l] += 0.5 * (log(SCHOOL[k].old_item_mat[i][j]) - olddelta[((i-1)*(i-2)/2+j)][l]) * (log(SCHOOL[k].old_item_mat[i][j]) - olddelta[((i-1)*(i-2)/2+j)][l]);
								else post_fb[l] += 0.5 * (log(0.0001) - olddelta[((i-1)*(i-2)/2+j)][l]) * (log(0.0001) - olddelta[((i-1)*(i-2)/2+j)][l]);
							}
					for(l = 0; l <= nCAT; l++) oldtau[((i-1)*(i-2)/2+j)][l] = 1.0 / Rgamma(post_fa[l], post_fb[l]);

					for(l = 0; l <= nCAT; l++){
						var_fix[l] = 1.0 / (1.0 / pr_var_delta + tvec[l] / oldtau[((i-1)*(i-2)/2+j)][l]);
						avg_fix[l] = 0.0;
						for(k = 1; k <= nSCHOOL; k++)
							if(t[k] == l){
								if(SCHOOL[k].old_item_mat[i][j] > 0.0001) avg_fix[l] += (1.0 / oldtau[((i-1)*(i-2)/2+j)][l]) * log(SCHOOL[k].old_item_mat[i][j]); 
								else avg_fix[l] += (1.0 / oldtau[((i-1)*(i-2)/2+j)][l]) * log(0.0001); 
							}
						avg_fix[l] *= var_fix[l];
						olddelta[((i-1)*(i-2)/2+j)][l] = avg_fix[l] + sqrt(var_fix[l]) * gasdev();
					}
				}

			if(iter > nburn && iter % thin == 0){
				gcount++;
				for(k = 1; k <= nSCHOOL; k++){
					sample_sigma[gcount][k] = sqrt(oldsigma[k]);
					fprintf(HUR, "%.4f ", sample_sigma[gcount][k]);
				}
				for(k = 0; k <= nCAT; k++)
					for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++){
						sample_tau[gcount][i][k] = sqrt(oldtau[i][k]);
						sample_delta[gcount][i][k] = olddelta[i][k];
						fprintf(JYW, "% .4f ", sample_delta[gcount][i][k]);
						fprintf(OUT, "%.4f ", sample_tau[gcount][i][k]);
					}
				for(k = 1; k <= nSCHOOL; k++)
					for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) sum_mu[i][k] += oldmu[i][k] / ((niter-nburn)/thin);
				for(i = 1; i <= nITEM; i++){
					for(k = 0; k <= nCAT; k++){
						sample_gamma[gcount][i][k] = oldgamma[i][k];
						sample_varphi[gcount][i][k] = sqrt(oldvarphi[i][k]);
						fprintf(JIN, "% .4f ", sample_gamma[gcount][i][k]);
						fprintf(PRT, "%.4f ", sample_varphi[gcount][i][k]);
					}
				}
				for(i = 1; i <= nSCHOOL; i++)
					for(j = 1; j <= nSCHOOL; j++) mu_dist[i][j] = 0.0;
				for(k = 1; k <= nITEM * (nITEM - 1) / 2; k++)
					for(i = 2; i <= nSCHOOL; i++)
						for(j = 1; j < i; j++) mu_dist[i][j] += (oldmu[k][i] - oldmu[k][j]) * (oldmu[k][i] - oldmu[k][j]);
				for(i = 2; i <= nSCHOOL; i++)
					for(j = 1; j < i; j++) mu_dist[j][i] = mu_dist[i][j];
				for(i = 1; i <= nSCHOOL; i++)
					for(j = 1; j <= nSCHOOL; j++) sum_mu_dist[i][j] += sqrt(mu_dist[i][j]) / ((niter-nburn)/thin);
				for(i = 2; i <= nSCHOOL; i++)
					for(j = 1; j < i; j++) fprintf(JJW, "%.4f ", sqrt(mu_dist[i][j]));

				fprintf(HUR, "\n"); fprintf(OUT, "\n"); fprintf(JYW, "\n");
				fprintf(JIN, "\n"); fprintf(PRT, "\n"); fprintf(JJW, "\n");
			}
		}
		fclose(HUR); fclose(JYW); fclose(OUT);
		fclose(JIN); fclose(PRT); fclose(JJW);

		frname[0] = 'R'; frname[1] = 'E'; frname[2] = 'S'; frname[3] = 'U'; frname[4] = 'L'; frname[5] = 'T'; frname[6] = '/';
		frname[7] = 's'; frname[8] = 'i'; frname[9] = 'm'; frname[12] = '_'; frname[14] = (char)(48+MM); frname[15] = '.';
		frname[16] = 'l'; frname[17] = 'o'; frname[18] = 'g'; frname[19] = '\0';
		for(a = 1; a <= nSCHOOL; a++){
			if(a < 10){frname[10] = (char)(48); frname[11] = (char)(a + 48);}
			else{frname[10] = (char)(a/10 + 48); frname[11] = (char)(a%10 + 48);}

			frname[13] = 'z'; JIN = fopen(frname, "a");
			frname[13] = 'b'; HUR = fopen(frname, "a");
			frname[13] = 't'; OUT = fopen(frname, "a");
			frname[13] = 'i'; JYW = fopen(frname, "a");
			frname[13] = 'h'; ASA = fopen(frname, "a");

			for(k = 1; k <= count[a]; k++){
				for(i = 1; i <= ncount[a]; i++)
					for(j = 1; j <= nDIM; j++) fprintf(JIN, "% .4f ", SCHOOL[a].sample_Zsamp[k][i][j]);
				fprintf(JIN, "\n");

				for(i = 1; i <= nITEM; i++)
					for(j = 1; j <= nDIM; j++) fprintf(JYW, "% .4f ", SCHOOL[a].sample_Zitem[k][i][j]);
				fprintf(JYW, "\n");

				for(i = 1; i <= nITEM; i++) fprintf(HUR, "% .4f ", SCHOOL[a].sample_beta[k][i]);
				fprintf(HUR, "\n");

				for(i = 1; i <= ncount[a]; i++) fprintf(OUT, "% .4f ", SCHOOL[a].sample_theta[k][i]);
				fprintf(OUT, "\n");

				fprintf(ASA, "%.4f\n", SCHOOL[a].sample_sigma[k]);
			}

			fclose(JIN); fclose(HUR); fclose(OUT); fclose(JYW); fclose(ASA);
		}

		// Calculate Mean and Variance of MCMC Estimators
		for(a = 1; a <= nSCHOOL; a++){
			for(i = 1; i <= count[a]; i++){
				SCHOOL[a].sum_sigma += SCHOOL[a].sample_sigma[i] / count[a];
				SCHOOL[a].var_sigma += SCHOOL[a].sample_sigma[i] * SCHOOL[a].sample_sigma[i] / (count[a] - 1);
				for(j = 1; j <= nITEM; j++){
					SCHOOL[a].sum_beta[j] += SCHOOL[a].sample_beta[i][j] / count[a];
					SCHOOL[a].var_beta[j] += SCHOOL[a].sample_beta[i][j] * SCHOOL[a].sample_beta[i][j] / (count[a] - 1);
				}
				for(j = 1; j <= ncount[a]; j++){
					SCHOOL[a].sum_theta[j] += SCHOOL[a].sample_theta[i][j] / count[a];
					SCHOOL[a].var_theta[j] += SCHOOL[a].sample_theta[i][j] * SCHOOL[a].sample_theta[i][j] / (count[a] - 1);
				}
				for(j = 1; j <= ncount[a]; j++)
					for(k = 1; k <= nDIM; k++){
						SCHOOL[a].sum_Zsamp[j][k] += SCHOOL[a].sample_Zsamp[i][j][k] / count[a];
						SCHOOL[a].var_Zsamp[j][k] += SCHOOL[a].sample_Zsamp[i][j][k] * SCHOOL[a].sample_Zsamp[i][j][k] / (count[a] - 1);
					}
				for(j = 1; j <= nITEM; j++)
					for(k = 1; k <= nDIM; k++){
						SCHOOL[a].sum_Zitem[j][k] += SCHOOL[a].sample_Zitem[i][j][k] / count[a];
						SCHOOL[a].var_Zitem[j][k] += SCHOOL[a].sample_Zitem[i][j][k] * SCHOOL[a].sample_Zitem[i][j][k] / (count[a] - 1);
					}
				for(j = 1; j <= nITEM; j++)
					for(k = 1; k <= nITEM; k++) SCHOOL[a].sample_item_mat[j][k] = 0.0;
				for(j = 2; j <= nITEM; j++)
					for(k = 1; k < j; k++)
						for(l = 1; l <= nDIM; l++) SCHOOL[a].sample_item_mat[j][k] += pow((SCHOOL[a].sample_Zitem[i][j][l] - SCHOOL[a].sample_Zitem[i][k][l]), 2.0);
				for(j = 2; j <= nITEM; j++)
					for(k = 1; k < j; k++) SCHOOL[a].sample_item_mat[k][j] = SCHOOL[a].sample_item_mat[j][k];
				for(j = 1; j <= nITEM; j++)
					for(k = 1; k <= nITEM; k++){
						SCHOOL[a].sum_item_mat[j][k] += SCHOOL[a].sample_item_mat[j][k] / count[a];
						SCHOOL[a].var_item_mat[j][k] += SCHOOL[a].sample_item_mat[j][k] * SCHOOL[a].sample_item_mat[j][k] / (count[a] - 1);
					}
			}
			SCHOOL[a].var_sigma -= SCHOOL[a].sum_sigma * SCHOOL[a].sum_sigma * count[a] / (count[a] - 1);
			for(i = 1; i <= nITEM; i++) SCHOOL[a].var_beta[i] -= SCHOOL[a].sum_beta[i] * SCHOOL[a].sum_beta[i] * count[a] / (count[a] - 1);
			for(i = 1; i <= ncount[a]; i++) SCHOOL[a].var_theta[i] -= SCHOOL[a].sum_theta[i] * SCHOOL[a].sum_theta[i] * count[a] / (count[a] - 1);
			for(i = 1; i <= ncount[a]; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[a].var_Zsamp[i][j] -= SCHOOL[a].sum_Zsamp[i][j] * SCHOOL[a].sum_Zsamp[i][j] * count[a] / (count[a] - 1);
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) SCHOOL[a].var_Zitem[i][j] -= SCHOOL[a].sum_Zitem[i][j] * SCHOOL[a].sum_Zitem[i][j] * count[a] / (count[a] - 1);
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nITEM; j++) SCHOOL[a].var_item_mat[i][j] -= SCHOOL[a].sum_item_mat[i][j] * SCHOOL[a].sum_item_mat[i][j] * count[a] / (count[a] - 1);
		}

		for(i = 1; i <= gcount; i++){
			for(k = 1; k <= nSCHOOL; k++){
				sum_sigma[k] += sample_sigma[i][k] / gcount;
				var_sigma[k] += sample_sigma[i][k] * sample_sigma[i][k] / (gcount - 1);
			}
			for(j = 1; j <= nITEM * (nITEM - 1) / 2; j++)
				for(k = 0; k <= nCAT; k++){
					sum_tau[j][k] += sample_tau[i][j][k] / gcount;
					sum_delta[j][k] += sample_delta[i][j][k] / gcount;
					var_tau[j][k] += sample_tau[i][j][k] * sample_tau[i][j][k] / (gcount - 1);
					var_delta[j][k] += sample_delta[i][j][k] * sample_delta[i][j][k] / (gcount - 1);
				}
			for(j = 1; j <= nITEM; j++)
				for(k = 0; k <= nCAT; k++){
					sum_gamma[j][k] += sample_gamma[i][j][k] / gcount;
					sum_varphi[j][k] += sample_varphi[i][j][k] / gcount;
					var_gamma[j][k] += sample_gamma[i][j][k] * sample_gamma[i][j][k] / (gcount - 1);
					var_varphi[j][k] += sample_varphi[i][j][k] * sample_varphi[i][j][k] / (gcount - 1);
				}
		}

		for(k = 0; k <= nSCHOOL; k++) var_sigma[k] -= sum_sigma[k] * sum_sigma[k] * gcount / (gcount - 1);
		for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++)
			for(k = 0; k <= nCAT; k++){
				var_tau[i][k] -= sum_tau[i][k] * sum_tau[i][k] * gcount / (gcount - 1);
				var_delta[i][k] -= sum_delta[i][k] * sum_delta[i][k] * gcount / (gcount - 1);
			}
		for(i = 1; i <= nITEM; i++)
			for(k = 0; k <= nCAT; k++){
				var_gamma[i][k] -= sum_gamma[i][k] * sum_gamma[i][k] * gcount / (gcount - 1);
				var_varphi[i][k] -= sum_varphi[i][k] * sum_varphi[i][k] * gcount / (gcount - 1);
			}

		// Save Parameter Estimates
		frname[0] = 'R'; frname[1] = 'E'; frname[2] = 'S'; frname[3] = 'U'; frname[4] = 'L'; frname[5] = 'T'; frname[6] = '/';
		frname[7] = 's'; frname[8] = 'u'; frname[9] = 'm'; frname[12] = '_'; frname[14] = (char)(48+MM); frname[15] = '.';
		frname[16] = 'l'; frname[17] = 'o'; frname[18] = 'g'; frname[19] = '\0';

		for(a = 1; a <= nSCHOOL; a++){
			if(a < 10){frname[10] = (char)(48); frname[11] = (char)(a + 48);}
			else{frname[10] = (char)(a/10 + 48); frname[11] = (char)(a%10 + 48);}

			frname[13] = 'z'; JIN = fopen(frname, "a");
			frname[13] = 'b'; HUR = fopen(frname, "a");
			frname[13] = 't'; OUT = fopen(frname, "a");
			frname[13] = 'i'; JYW = fopen(frname, "a");
			frname[13] = 'd'; PRT = fopen(frname, "a");

			for(i = 1; i <= nITEM; i++) fprintf(HUR, "%.4f ", SCHOOL[a].sum_beta[i]); fprintf(HUR, "\n");
			for(i = 1; i <= nITEM; i++) fprintf(HUR, "%.4f ", SCHOOL[a].var_beta[i]); fprintf(HUR, "\n");
			for(i = 1; i <= nITEM; i++) fprintf(HUR, "%.4f ", SCHOOL[a].acc_beta[i]); fprintf(HUR, "\n");

			for(i = 1; i <= ncount[a]; i++) fprintf(OUT, "%.4f ", SCHOOL[a].sum_theta[i]); fprintf(OUT, "\n");
			for(i = 1; i <= ncount[a]; i++) fprintf(OUT, "%.4f ", SCHOOL[a].var_theta[i]); fprintf(OUT, "\n");
			for(i = 1; i <= ncount[a]; i++) fprintf(OUT, "%.4f ", SCHOOL[a].acc_theta[i]); fprintf(OUT, "\n");

			for(i = 1; i <= ncount[a]; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JIN, "%.4f ", SCHOOL[a].sum_Zsamp[i][j]); fprintf(JIN, "\n");
			for(i = 1; i <= ncount[a]; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JIN, "%.4f ", SCHOOL[a].var_Zsamp[i][j]); fprintf(JIN, "\n");
			for(i = 1; i <= ncount[a]; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JIN, "%.4f ", SCHOOL[a].acc_Zsamp[i]); fprintf(JIN, "\n");

			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JYW, "%.4f ", SCHOOL[a].sum_Zitem[i][j]); fprintf(JYW, "\n");
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JYW, "%.4f ", SCHOOL[a].var_Zitem[i][j]); fprintf(JYW, "\n");
			for(i = 1; i <= nITEM; i++)
				for(j = 1; j <= nDIM; j++) fprintf(JYW, "%.4f ", SCHOOL[a].acc_Zitem[i]); fprintf(JYW, "\n");

			for(i = 2; i <= nITEM; i++)
				for(j = 1; j < i; j++) fprintf(PRT, "%.4f ", SCHOOL[a].sum_item_mat[i][j]); fprintf(PRT, "\n");
			for(i = 2; i <= nITEM; i++)
				for(j = 1; j < i; j++) fprintf(PRT, "%.4f ", SCHOOL[a].var_item_mat[i][j]); fprintf(PRT, "\n");

			fclose(JIN); fclose(HUR); fclose(OUT); fclose(JYW); fclose(PRT);
		}

		frname[0] = 'R'; frname[1] = 'E'; frname[2] = 'S'; frname[3] = 'U'; frname[4] = 'L'; frname[5] = 'T'; frname[6] = '/';
		frname[7] = 's'; frname[8] = 'u'; frname[9] = 'm'; frname[10] = '_'; frname[12] = (char)(48+MM); frname[13] = '.';
		frname[14] = 'l'; frname[15] = 'o'; frname[16] = 'g'; frname[17] = '\0';

		frname[11] = 'm'; JIN = fopen(frname, "a");
		frname[11] = 's'; HUR = fopen(frname, "a");
		frname[11] = 'l'; JYW = fopen(frname, "a");
		frname[11] = 'u'; OUT = fopen(frname, "a");
		frname[11] = 'g'; ASA = fopen(frname, "a");
		frname[11] = 'p'; PRT = fopen(frname, "a");
		frname[11] = 'h'; IMS = fopen(frname, "a");
		frname[11] = 'a'; JJW = fopen(frname, "a");

		for(i = 1; i <= nSCHOOL; i++) fprintf(HUR, "%.4f ", sum_sigma[i]); fprintf(HUR, "\n");
		for(i = 1; i <= nSCHOOL; i++) fprintf(HUR, "%.4f ", var_sigma[i]); fprintf(HUR, "\n");

		for(k = 0; k <= nCAT; k++){
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) fprintf(OUT, "%.4f ", sum_tau[i][k]); fprintf(OUT, "\n");
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) fprintf(OUT, "%.4f ", var_tau[i][k]); fprintf(OUT, "\n");

			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) fprintf(JYW, "% .4f ", sum_delta[i][k]); fprintf(JYW, "\n");
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) fprintf(JYW, "% .4f ", var_delta[i][k]); fprintf(JYW, "\n");
		}
		for(k = 1; k <= nSCHOOL; k++){
			for(i = 1; i <= nITEM * (nITEM - 1) / 2; i++) fprintf(JIN, "% .4f ", sum_mu[i][k]); fprintf(JIN, "\n");
		}
		for(k = 0; k <= nCAT; k++){
			for(i = 1; i <= nITEM; i++) fprintf(ASA, "% .4f ", sum_gamma[i][k]); fprintf(ASA, "\n");
			for(i = 1; i <= nITEM; i++) fprintf(ASA, "% .4f ", var_gamma[i][k]); fprintf(ASA, "\n");
			for(i = 1; i <= nITEM; i++) fprintf(PRT, "%.4f ", sum_varphi[i][k]); fprintf(PRT, "\n");
			for(i = 1; i <= nITEM; i++) fprintf(PRT, "%.4f ", var_varphi[i][k]); fprintf(PRT, "\n");
		}

		for(k = 1; k <= nSCHOOL; k++) fprintf(IMS, "%.4f ", SCHOOL[k].sum_sigma); fprintf(IMS, "\n");
		for(k = 1; k <= nSCHOOL; k++) fprintf(IMS, "%.4f ", SCHOOL[k].var_sigma); fprintf(IMS, "\n");
		for(i = 1; i <= nSCHOOL; i++){
			for(j = 1; j <= nSCHOOL; j++) fprintf(JJW, "%.4f ", sum_mu_dist[i][j]); fprintf(JJW, "\n");
		}
		
		fclose(JIN); fclose(HUR); fclose(JYW); fclose(OUT);
		fclose(ASA); fclose(PRT); fclose(IMS); fclose(JJW);
	}

	/*
	free_ivector(ncount, 1, nSCHOOL);
	free_dvector(jump_Z, 0, nITEM);

	for(k = 0; k <= nSCHOOL; k++){
		for(i = 0; i <= nMAX; i++) free(SCHOOL[k].dataset[i]);
		for(i = 0; i <= nITEM; i++){
			for(a = 0; a <= nMAX; a++) free(SCHOOL[k].Y[i][a]);
			free(SCHOOL[k].Y[i]);
		}
 		for(i = 0; i <= nMAX; i++){
			for(a = 0; a <= nITEM; a++) free(SCHOOL[k].U[i][a]);
			free(SCHOOL[k].U[i]);
		}
		for(i = 0; i <= nMAX; i++){free(SCHOOL[k].old_Zsamp[i]); free(SCHOOL[k].new_Zsamp[i]);}
		for(i = 0; i <= nITEM; i++){free(SCHOOL[k].old_Zitem[i]); free(SCHOOL[k].new_Zitem[i]);}
 		for(i = 0; i <= (niter-nburn)/thin; i++){
			for(j = 0; j <= nMAX; j++)  free(SCHOOL[k].sample_Zsamp[i][j]);
			for(j = 0; j <= nITEM; j++) free(SCHOOL[k].sample_Zitem[i][j]);
			free(SCHOOL[k].sample_beta[i]); free(SCHOOL[k].sample_theta[i]);
			free(SCHOOL[k].sample_Zsamp[i]); free(SCHOOL[k].sample_Zitem[i]);
		}
		for(i = 0; i <= nMAX; i++){free(SCHOOL[k].sum_Zsamp[i]); free(SCHOOL[k].var_Zsamp[i]);}
		for(i = 0; i <= nITEM; i++){free(SCHOOL[k].sum_Zitem[i]); free(SCHOOL[k].var_Zitem[i]);}
		for(i = 0; i <= nITEM; i++){free(SCHOOL[k].sum_item_mat[i]); free(SCHOOL[k].var_item_mat[i]);}
		for(i = 0; i <= nITEM; i++){free(SCHOOL[k].old_item_mat[i]); free(SCHOOL[k].new_item_mat[i]); free(SCHOOL[k].sample_item_mat[i]);}
		free(SCHOOL[k].old_item_mat); free(SCHOOL[k].new_item_mat);
		free(SCHOOL[k].oldbeta); free(SCHOOL[k].newbeta);
		free(SCHOOL[k].oldtheta); free(SCHOOL[k].newtheta);
		free(SCHOOL[k].count_item); free(SCHOOL[k].count_samp);
		free(SCHOOL[k].Y); free(SCHOOL[k].U); free(SCHOOL[k].dataset);
		free(SCHOOL[k].old_Zsamp); free(SCHOOL[k].new_Zsamp);
		free(SCHOOL[k].old_Zitem); free(SCHOOL[k].new_Zitem);
		free(SCHOOL[k].sample_beta); free(SCHOOL[k].sample_theta);
		free(SCHOOL[k].sum_beta); free(SCHOOL[k].var_beta); free(SCHOOL[k].acc_beta);
		free(SCHOOL[k].sum_theta); free(SCHOOL[k].var_theta); free(SCHOOL[k].acc_theta);
		free(SCHOOL[k].sample_Zsamp); free(SCHOOL[k].sample_Zitem); free(SCHOOL[k].sample_item_mat);
		free(SCHOOL[k].sum_Zsamp); free(SCHOOL[k].var_Zsamp); free(SCHOOL[k].acc_Zsamp);
		free(SCHOOL[k].sum_Zitem); free(SCHOOL[k].var_Zitem);
		free(SCHOOL[k].sum_item_mat); free(SCHOOL[k].var_item_mat);
		free(SCHOOL[k].sample_sigma); free(SCHOOL[k].mean_Z);
	}
	free(SCHOOL);

	for(i = 0; i <= (niter-nburn)/thin; i++){
		for(j = 0; j <= nITEM * nDIM; j++){
			free(sample_tau[i][j]);
			free(sample_delta[i][j]);
		}
		for(j = 0; j <= nITEM; j++){
			free(sample_gamma[i][j]);
			free(sample_varphi[i][j]);
		}
		free(sample_tau[i]);
		free(sample_delta[i]);
		free(sample_gamma[i]);
		free(sample_varphi[i]);
	}
	free(sample_tau);
	free(sample_delta);
	free(sample_gamma);
	free(sample_varphi);
	free_dmatrix(sample_sigma, 1, (niter-nburn)/thin, 1, nSCHOOL);
	free_dmatrix(sum_mu, 1, nITEM * nDIM, 0, nSCHOOL);
	free_dmatrix(sum_tau, 1, nITEM * nDIM, 0, nCAT);
	free_dmatrix(var_tau, 1, nITEM * nDIM, 0, nCAT);
	free_dvector(sum_sigma, 1, nSCHOOL);
	free_dvector(var_sigma, 1, nSCHOOL);
	free_dmatrix(sum_delta, 1, nITEM * nDIM, 0, nCAT);
	free_dmatrix(var_delta, 1, nITEM * nDIM, 0, nCAT);
	free_dmatrix(sum_gamma, 1, nITEM, 0, nCAT);
	free_dmatrix(var_gamma, 1, nITEM, 0, nCAT);
	free_dmatrix(sum_varphi, 1, nITEM, 0, nCAT);
	free_dmatrix(var_varphi, 1, nITEM, 0, nCAT);

	free_dvector(post_ra, 1, nSCHOOL); free_dvector(post_rb, 1, nSCHOOL);
	free_dvector(post_fa, 0, nCAT); free_dvector(post_fb, 0, nCAT);
	free_dvector(school_a, 0, nCAT); free_dvector(school_b, 0, nCAT);

	free_ivector(t, 1, nSCHOOL);
	free_ivector(tvec, 0, nCAT);
	free_dvector(oldsigma, 1, nSCHOOL);
	free_dmatrix(olddelta, 1, nITEM * nDIM, 0, nCAT);
	free_dmatrix(oldtau, 1, nITEM * nDIM, 0, nCAT);
	free_dmatrix(oldmu, 1, nITEM * nDIM, 0, nSCHOOL);
	free_dmatrix(oldgamma, 1, nITEM, 0, nCAT);
	free_dmatrix(oldvarphi, 1, nITEM, 0, nCAT);

	free_dvector(var_fix, 0, nCAT);
	free_dvector(avg_fix, 0, nCAT);
	free_dvector(var_ran, 1, nSCHOOL);
	free_dvector(avg_ran, 1, nSCHOOL);
	free_dvector(avg_beta, 0, nCAT);
	free_dvector(var_beta, 0, nCAT);

	free_dmatrix(mu_dist, 1, nSCHOOL, 1, nSCHOOL);
	free_dmatrix(sum_mu_dist, 1, nSCHOOL, 1, nSCHOOL);

	free_ivector(count, 1, nSCHOOL);

	free_dvector(sample_samp_like, 1, nMAX);
	free_dvector(new_samp_distance, 1, nMAX);
	free_dvector(old_samp_distance, 1, nMAX);

	free_dvector(sample_item_like, 1, nMAX);
	free_dmatrix(new_item_distance, 1, nITEM, 1, nITEM);
	free_dmatrix(old_item_distance, 1, nITEM, 1, nITEM);
	*/

	return 0;
}
