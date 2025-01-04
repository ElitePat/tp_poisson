/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

#define MAX(i, j) (((i) > (j)) ? (i) : (j)) // très utile !

// Stockage GB (General Band) en priorite colonne pour matrice Poisson 1D
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
	/* AB represente la matrice (à remplir)
	lab: nombre de diagonales de la matrice (diagonale)
	la: taille maximale des diagonales
	kv: offsets de securité (nombre de lignes supplementaires)

	Soit pour une matrice de Poisson 1D on aurait l'appel:
	set_GB_operator_colMajor_poisson1D(AB, 3, la, 1);
	*/

	// Variables 
	int size = (*lab) * (*la) + ((*la) * (*kv)); // nombre total d'elemnts AB
	int ld = (*lab) + (*kv); // nombre de lignes de AB
	int cntr = (ld / 2) + 1; // la ligne centrale (à remplir avec des 2) 
	// normalement cntr = 2 puisque matrice Poisson 1D
	int i;

	// nous allons passer par 3 boucles:
	for(i=0;i<size;i++){ // on met tout à zéro
		AB[i] = 0;
	}
	for(i=cntr;i<size;i+ld){ // on rempli la ligne centrale
		AB[i] = 2;
	}
	for(int j=1;j<ld;j+2){
		for(i=j;i<size-1;i+ld){ // size-1 pour ne pas faire la dernière valeur
			if(i != 1){ // on saute aussi l'elemtent a(1,2) de la
				AB[i] = -1;
			}
		}
	}
}

// Identité (peut servir pour les méthodes de validation)
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
	// on applique le même raisonement que pour la fonction precedente:

	// Variables 
	int size = (*lab) * (*la) + ((*la) * (*kv)); // nombre total d'elemnts AB
	int ld = (*lab) + (*kv); // nombre de lignes de AB (normalement 2)
	int cntr = ld - 1; // car matrice identité !
	int i;

	for(i=0;i<size;i++){ // le reste de la matrice
		AB[i] = 0;
	}
	for(i=cntr;i<size;i+ld){ // on saute de colone en colone
		AB[i] = 1;// la ligne principale (à remplir avec des 1)
	}
}


void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
	// correction
	int jj;
	RHS[0]= *BC0;
	RHS[(*la)-1]= *BC1;
	for (jj=1;jj<(*la)-1;jj++){
		RHS[jj]=0.0;
	}
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
	// correction
	int jj;
	double h, DELTA_T;
	DELTA_T=(*BC1)-(*BC0);
	for (jj=0;jj<(*la);jj++){
		EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
	}
}  


void set_grid_points_1D(double* x, int* la){
	// correction
	int jj;
	double h;
	h=1.0/(1.0*((*la)+1));
	for (jj=0;jj<(*la);jj++){
		x[jj]=(jj+1)*h;
	}
}

double relative_forward_error(double* x, double* y, int* la){
	/*la : taille de la matrice 
	y : matrice resultat (On suppose que x et y sont carrées)
	On a pas assez d'info pour faire appel aux fonctions lapack qui 
	calculent la norme matricielle.
	Alors on implemente "à la main" la definition de la norme infinie. */

	// x - y:
	double r[(*la)];
	for(int a=0;a<(*la)*(*la);a++){
		r[a] = x[a] - y[a];
	}

	// Norme infinie de (x-y) et de x
	double xmax, xsome, rmax, rsome;
	// on calcule les 2 en même temps (traitement identique)
	for(int i=0;i<(*la);i++){
		xsome = 0;
		rsome = 0;
		for(int j=0;j<(*la);j++){
			// à priori fabs() retourne un double
			xsome += fabs(x[i]);
			rsome += fabs(y[i]);
		}
		xmax = MAX(xmax,xsome);
		rmax = MAX(rmax,rsome);
	}

	//resultat
	return rmax / xmax;
}

int indexABCol(int i, int j, int *lab){
	return j*(*lab)+i;
}

// factorisatoin LU pour matrices tridiagonales
int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
	/*
	la: longueur diagonale
	n: 
	kl: matrice triangulaire inferieure
	ku: matrice triangulaire superieure
	lab: à priori ègal a 3 car tridigonale
	AB: la matrice de départ (en GB)
	*/	
	return *info;
}
