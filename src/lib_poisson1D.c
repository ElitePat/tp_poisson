/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

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
	int size = lab * la + (la * kv); // nombre total d'elemnts AB
	int ld = lab + kv; // nombre de lignes de AB
	int cntr = lab - kv; // la ligne centrale (à remplir avec des 2) 
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

// Identité 
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
	// on applique le même raisonement que pour la fonction precedente:

	// Variables 
	int size = lab * la + (la * kv); // nombre total d'elemnts AB
	int ld = lab + kv; // nombre de lignes de AB
	int cntr = (int)((ld / 2.0) + 1); // la ligne centrale (à remplir avec des 1)
	int i;

	/* en vrai de vrai, à voir ================================================= */
}

// 
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

// la = |vevteur|
double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
	return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
	return *info;
}
