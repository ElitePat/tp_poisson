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
	return 4 * pow(sin(M_PI * (1.0/((*la) + 1) / 2.0)),2);
}

double richardson_alpha_opt(int *la){
	return 2.0 / (eigmin_poisson1D(&la) + eigmax_poisson1D(&la));
}
/*	On notera que l'alpha optimal pour Jaccobi et Gaus-Siedel est de 1/2
*/

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	/* AB : matrice poisson1D 
	RHS : matrice dense (à priori)
	X : matrice/vecteur solution
	tol : la tolerance
	resvec: vecteur contenant le résidu */

	double *Y = malloc((*la)*sizeof(double));

	cblas_dcopy((*la),Y,RHS); // arguments like N, X, incX, Y, incY
	double nomrY = cblas_nomr2(Y);
	cblas_dgbmv(CblasColMajor,CblasNoTrans,(*la),(*la),(*kl),(*ku),-1.0,AB,(*lab),X,1,1.0,Y,1);
	double nomr_res = cblas_dnorm(Y);
	double res = nomr_res / nomrY;

	while((res > (*tol)) && ((*nbite) < (*maxit))){
		cblas_daxpy((*la),(*alpha_rich),Y,1,X,1);
		resvec[(*nbite)] = res;
		(*nbite)++;
	}

	free(Y);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
	
}

// le truc general
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

}

