/******************************************/
/* tp_env.c                               */
/* This file contains a main function to  */
/* test the environment of compilation    */
/******************************************/
#include "tp_env.h"

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
	printf("--------- Test environment of execution for Practical exercises of Numerical Algorithmics ---------\n\n");

	/*
	printf("Test de la méthode dgbmv de BLAS\n");
	// declaration des varaibles
	double alpha=1 ,beta=2; // scalaires
	double *GB, *ID, *x, *y; // matrice poisson 1D dans GB, matrice Id en ID et vecteurs
	int un=1, trois=3, N=4; // N=taille matrices (pour l'instant carrées)


	// phase d'allocation mémoire
	x = (double *) malloc(sizeof(double)*N);
	y = (double *) malloc(sizeof(double)*N);
	GB = (double *) malloc(sizeof(double)*N*N);

	// instantiation des variables
	set_GB_operator_colMajor_poisson1D_Id(ID,&un,&N,&un);
	set_GB_operator_colMajor_poisson1D(GB, &trois, &N, &un);
	// afichage de la matrice
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			printf(" %lf",GB[N*i+j]);
		}
		printf("\n");
	}

	for(int i=0;i<N;i++){
		x[i] = 0;
		y[i] = 0;
	} // pour l'instant on fait simple


	// appel a cblas_dgbmv qui faira le calcul suivant: y := alpha*GB*x + beta*y
	cblas_dgbmv(CblasColMajor,CblasNoTrans,N,N,1,1,alpha,GB,N,x,N,beta,y,N);
	// c'est le vecteur y qui contient le resultat de l'appel

	// affichage du resultat y:
	for(int k=0;k<N;k++){
		printf("| %lf |\n",y[k]);
	}

	free(x);
	free(y);
	free(GB);
	*/

	
	printf("The exponantial value is e = %f \n",M_E);
	printf("The maximum single precision value from values.h is maxfloat = %e \n",MAXFLOAT);
	printf("The maximum single precision value from float.h is flt_max = %e \n",FLT_MAX);
	printf("The maximum double precision value from float.h is dbl_max = %e \n",DBL_MAX);
	printf("The epsilon in single precision value from float.h is flt_epsilon = %e \n",FLT_EPSILON);
	printf("The epsilon in double precision value from float.h is dbl_epsilon = %e \n",DBL_EPSILON);

	printf("\n\n Test of ATLAS (BLAS/LAPACK) environment \n");

	double x[5], y[5];
	int ii;
	for (ii=0;ii<5;ii++){
		x[ii]=ii+1; y[ii]=ii+6;
		printf("x[%d] = %lf, y[%d] = %lf\n",ii,x[ii],ii,y[ii]);
	}

	printf("\nTest DCOPY y <- x \n");
	cblas_dcopy(5,x,1,y,1);
	for (ii=0;ii<5;ii++){
		printf("y[%d] = %lf\n",ii,y[ii]);
	}

	printf("\n\n--------- End -----------\n");
}
