/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// Stockage GB en priorite colonne pour la matrice de Poisson 1D
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
	// AB est vide !
	/*  Pour une matrice de Poisson 1D on aurait
	set_GB_operator_colMajor_poisson1D(AB, 3, la, 1);
	*/
	int size = (*lab) * (*la) + (*kv); // nombre total d'elemnts AB
	int ld = (*lab) + (*kv); // nombre de lignes de AB
	int cntr = (*lab) - (*kv);
	int kkv = (*kv); // on prend la première valeur
	for(int i=0;i<size;i++){
		if(i%ld == 0){
			AB[i] = 0;
		}else if((i%size = kkv)%){}
	}
}

// Identité
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
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
