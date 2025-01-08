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
	return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
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

	cblas_dcopy((*la),Y,0,RHS,0); // arguments like N, X, incX, Y, incY
	double nomrY = cblas_dnrm2((*la),Y,0);
	cblas_dgbmv(CblasColMajor,CblasNoTrans,(*la),(*la),(*kl),(*ku),-1.0,AB,(*lab),X,1,1.0,Y,1); // alpha*A*x + beta*y
	double nomr_res = cblas_dnrm2((*la),Y,0);
	double res = nomr_res / nomrY;

	while((res > (*tol)) && ((*nbite) < (*maxit))){
		cblas_daxpy((*la),(*alpha_rich),Y,1,X,1); // Ax + y
		resvec[(*nbite)] = res;
		(*nbite)++;
	}

	free(Y);
}

/* Dans les deux focntions qui suivent nous n'allos pas inverser MB, cela sera fait plus tard!*/

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
	/*
	Jaccobi utilise la diagonale de A [selon la formule x(k+1) = x(k) + M^(-1) * (b-A*x(k))]
	la: taille de la diagonale
	lab: nombre de lignes
	AB: matrice poisson 1D (?)
	MB: matrice diagonale de taille lab*la (en GB) => diag(A)
	ku: nombre de diagonales superieures
	kl: nombre de diagonales inferieures
	kv: taille de l'offset
	*/

	for(int i=0;i<(*la);i++){ // pour chaque colonne
		for(int j=0;j<((*lab)+(*kv));j++){ // pour chaque ligne
			if(j==((*ku)+(*kv))){ // on cible la diagonale
				MB[(*la)*i+j] = AB[(*la)*i+j]; // diag(A)
			}else{
				MB[(*la)*i+j] = 0; // pour le reste c'est zéro
			}
		}
	}

}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
	/*
	la: taille de la diagonale
	lab: nombre de lignes
	AB: matrice poisson 1D (?)
	MB: matrice diagonale de taille lab*la (en GB) => diag(A)-E
	ku: nombre de diagonales superieures
	kl: nombre de diagonales inferieures
	kv: taille de l'offset
	*/
	for(int i=0;i<(*la);i++){ // pour chaque colonne
		for(int j=0;j<((*lab)+(*kv));j++){ // pour chaque ligne
			if(j==((*ku)+(*kv))){ // on cible la diagonale et la diagonale inferieure
				MB[(*la)*i+j] = AB[(*la)*i+j]; // diag(A)
			}else if (j==((*ku)+(*kv)+1)){
				MB[(*la)*i+j] = (-1) * AB[(*la)*i+j]; // -E
			}else{
				MB[(*la)*i+j] = 0; // pour le reste c'est zéro
			}
		}
	}
}

// la definition de la méthode de Richardson: [x^(k+1) = x^k + M^(-1) * (b - A * x^k)]
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
	/*
	lab: nombre de lignes de AB
	la: taille de la diagonale de AB
	les arguments sont les mêmes que pour la fonction richardson_alpha() sauf
	MB: matrice d'itération (en GB) -> defini la méthode itérative de Richardson qu'on utilise
	RHS: contient le residu (Right Hand Side)
	*/

	double *B = malloc((*la)*sizeof(double));
	cblas_dcopy((*la),RHS,0,B,0); // on copie RHS dans B
	// on va utiliser RHS et laisser B constant

	// r^0 = (b - A * x^0) le residu 
	cblas_dgbmv(CblasColMajor,CblasNoTrans,(*la),(*lab),(*kl),(*ku),1.0,AB,(*la),X,1,-1.0,RHS,1); // RHS = (1.0 * AB * X + -1.0 * RHS)

	// Condition d'arret: (|| r^(k+1) || < eps) ou nb d'iterations maxit
	double norm_res;
	while((norm_res > (*tol)) || ((*nbite) < (*maxit))){
		/* Nous allons pas inverser la matrice mais nous allons faire ceci:
		x^(k+1) = x^k + M^(-1) * (b - A * x^k)
			<=> M*x^(k+1) = M*x^k + (b-A*x^k)
			<=> MB*X = MB*X + (AB*X + RHS=B)*/
		cblas_dcopy((*la),B,0,RHS,0); // on copie B dans RHS
		cblas_dgbmv(CblasColMajor,CblasNoTrans,(*la),(*lab),(*kl),(*ku),1.0,AB,(*la),X,1,-1.0,RHS,1); // RHS = (1.0 * AB * X + -1.0 * RHS)
		norm_res = cblas_dnrm2((*la),RHS,0); // calcul de la norme de RHS^(k+1)
		cblas_dgbmv(CblasColMajor,CblasNoTrans,(*la),(*lab),(*kl),(*ku),1.0,MB,(*la),X,1,-1.0,RHS,1); // RHS = (1.0 * MB * X + 1.0 * RHS)
		/* Au lieu de faire dgbsv() ici on aurait pu resoudre l'equation avec dgbtrftridiag()
		 à cause de MB qui est tridiadonale */
		//dgbsv (fact LU et LU*x=y)
		int info, *ipiv;
		dgbsv_((*la),(*kl),(*ku),1,MB,(*la),ipiv,RHS,(*la),&info);
		if(info!=0){printf("Attention DGBSV: info = %d",info);}
		cblas_dcopy((*la),RHS,0,X,0); // on copie RHS dans X
		(*nbite)++;
	}

	free(B);
}

