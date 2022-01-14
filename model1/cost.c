double cost_Zsamp(double **Z, int sh){

	int i, j, k, l;
	double log_like = 0.0, dist;

	#pragma omp parallel for private(i, j, k, l, dist) default(shared)
	for(k = 1; k <= nITEM; k++){
		for(i = 2; i <= ncount[sh]; i++)
			for(j = 1; j < i; j++){
				for(dist = 0.0, l = 1; l <= nDIM; l++) dist += pow((Z[i][l] - Z[j][l]), 2.0);
				dist = sqrt(dist);
				if(SCHOOL[sh].Y[k][i][j] == 1) log_like += -log(1.0 + exp(-(SCHOOL[sh].oldbeta[k] - dist)));
				else log_like += -log(1.0 + exp(SCHOOL[sh].oldbeta[k] - dist));
			}
	}

	return log_like;
}

double cost_Zitem(double **Z, int sh){

	int i, j, k, l;
	double log_like = 0.0, dist;

	#pragma omp parallel for private(i, j, k, l, dist) default(shared)
	for(k = 1; k <= ncount[sh]; k++){
		for(i = 2; i <= nITEM; i++)
			for(j = 1; j < i; j++){
				for(dist = 0.0, l = 1; l <= nDIM; l++) dist += pow((Z[i][l] - Z[j][l]), 2.0);
				dist = sqrt(dist);
				if(SCHOOL[sh].U[k][i][j] == 1) log_like += -log(1.0 + exp(-(SCHOOL[sh].oldtheta[k] - dist)));
				else log_like += -log(1.0 + exp(SCHOOL[sh].oldtheta[k] - dist));
			}
	}

	return log_like;
}

double cost_beta(int item, double beta, int sh){

	int i, j, l;
	double log_like = 0.0, dist;

	//#pragma omp parallel for private(i, j, l, dist) default(shared)
	for(i = 2; i <= ncount[sh]; i++)
		for(j = 1; j < i; j++){
			for(dist = 0.0, l = 1; l <= nDIM; l++) dist += pow((SCHOOL[sh].old_Zsamp[i][l] - SCHOOL[sh].old_Zsamp[j][l]), 2.0);
			dist = sqrt(dist);
			if(SCHOOL[sh].Y[item][i][j] == 1) log_like += -log(1 + exp(-(beta - dist)));
			else log_like += -log(1 + exp(beta - dist));
		}

	return log_like;
}

double cost_theta(int sample, double theta, int sh){

	int i, j, l;
	double log_like = 0.0, dist;

	//#pragma omp parallel for private(i, j, l, dist) default(shared)
	for(i = 2; i <= nITEM; i++)
		for(j = 1; j < i; j++){
			for(dist = 0.0, l = 1; l <= nDIM; l++) dist += pow((SCHOOL[sh].old_Zitem[i][l] - SCHOOL[sh].old_Zitem[j][l]), 2.0);
			dist = sqrt(dist);
			if(SCHOOL[sh].U[sample][i][j] == 1) log_like += -log(1 + exp(-(theta - dist)));
			else log_like += -log(1 + exp(theta - dist));
		}

	return log_like;
}

