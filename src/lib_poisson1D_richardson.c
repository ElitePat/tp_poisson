/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

#define MAX(i, j) (((i) > (j)) ? (i) : (j)) // très utile !

// remplir le tableau des valeurs propres
void eig_poisson1D(double* eigval, int *la){
	double h = (1.0 / ((double)(*la) + 1.0));
	double sinT;
	for(int k=1;k<=(*la);k++){
		sinT = sin(((double)k * M_PI * h)/ 2.0);
		eigval[k-1] = 4 * sinT * sinT;
	}
	/* On met k-1 car on compte les valeurs de 1 à K
	et en C on commence à 0 */
}

double eigmax_poisson1D(int *la){
	// on peut reutiliser la fonction précedente
	// où utiliser la formule directe
	return 4 * pow(sin( (*la) * M_PI * (1.0/((*la) + 1) / 2.0)),2);
}

double eigmin_poisson1D(int *la){
	return 0;
}

double richardson_alpha_opt(int *la){
	return 0;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

