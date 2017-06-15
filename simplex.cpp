#include"simplex.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

static const double eps = 1e-4;

typedef struct tableau {
	int m, n;
	double A[MAX][MAX];
} Tableau;

int igual(double x, double y) {
	return fabs(x - y) < eps;
}

void mostra_A(double A[][MAX], int m, int n) {
	int i, j;
	printf("\nMatriz A: \n");
	for (i = 1; i <= m; ++i) {
		for (j = 1; j <= n; ++j) {
			printf("%lf\t", A[i][j]);
		}
		printf("\n");
	}
}

void mostra_z(double z[], int n) {
	int i;
	printf("\nCoeficientes de z: \n");
	for (i = 1; i <= n; ++i) {
		printf("%lf ", z[i]);
	}
	printf("\n");
}

void mostra_b(double b[], int m) {
	int i, j;
	printf("\nConteudo de b: \n");
	for (i = 1; i <= m; ++i) {
		printf("%lf ", b[i]);
	}
	printf("\n");
}

int base_viavel(double b[], int m) {
	int i, j;

	for (i = 1; i <= m; ++i) {
		if (b[i] < 0)
			return 0;
	}

	return 1;
}

Tableau *cria_tableau(int m, int n, double z[], double A[][MAX], double b[]) {
	int i, j;
	Tableau *tab = (Tableau*) malloc(sizeof(Tableau));
	tab->m = m + 1;
	tab->n = n + 1;

	for (i = 0; i <= m; ++i)
		for (j = 0; j <= n; ++j)
			tab->A[i][j] = A[i][j];

	for (i = 1; i <= m; ++i)
		tab->A[i][0] = b[i];
	for (j = 1; j <= n; ++j)
		tab->A[0][j] = -z[j];

	return tab;
}

void mostra_tableau(Tableau *tab) {
	int i, j;
	printf("\n");
	for (i = 0; i < tab->m; ++i) {
		for (j = 0; j < tab->n; ++j) {
			printf("%.1lf\t", tab->A[i][j]);
		}
		printf("\n");
	}
}
//Adiciona variaveis de folga em cada restrição
void add_folgas(Tableau *tab) {
	int i, j, n = tab->n;
	tab->n += tab->m - 1;
	for (i = 1; i < tab->m; ++i) {
		for (j = n; j < tab->n; ++j) {
			if ((j - n) == i - 1)
				tab->A[i][j] = 1;
			else
				tab->A[i][j] = 0;
		}
	}
}

int coluna_pivo(Tableau *tab) {
	int i;
	for (i = 1; i < tab->n; ++i)
		if (tab->A[0][i] + eps < 0)
			return i;

	return -1;
}

//Encontra a linha pivô: menor razão p_i/x_v_i.
int linha_pivo(Tableau *tab, int col_pivo) {
	double razao;
	int i, lin_pivo;
	double min = INFINITY;

	for (i = 1; i < tab->m; ++i) {
		if (tab->A[i][col_pivo] - eps > 0) {
			razao = tab->A[i][0] / tab->A[i][col_pivo];
			printf("\nRazão: %lf", razao);
			if (min >= razao) { //Em caso de empate, permanece o de menor índice
				min = razao;
				printf("\nMin: %lf", min);
				lin_pivo = i;
			}
		}
	}
	//Se não tem elemento positivo na coluna pivô
	if (min == INFINITY)
		return -1;

	return lin_pivo;
}

//Usa regra de Bland para escolher a variavel que entra e a que sai
void pivoteamento(Tableau *tab, int lin_pivo, int col_pivo) {
	int i, j;
	double pivo;
	pivo = tab->A[lin_pivo][col_pivo];
	printf("\nPivô: %lf\n", pivo);

	//Nova linha pivô
	for (j = 0; j < tab->n; ++j) {
		tab->A[lin_pivo][j] = tab->A[lin_pivo][j] / pivo;
	}
	// Calcula as novas linhas
	for (i = 0; i < tab->m; ++i) {
		float tmp = tab->A[i][col_pivo];
		if (i != lin_pivo) {
			for (j = 0; j < tab->n; ++j) {
				tab->A[i][j] += tab->A[lin_pivo][j] * -tmp;
			}
		}
	}
	mostra_tableau(tab);
}
void vetor_otimo(Tableau *tab, double x[]) {
	int i, j, k = 1;
	for (j = 1; j < tab->n; ++j) {
		if (!igual(tab->A[0][j], 0)) {
			x[k++] = 0;
		} else {
			for (i = 1; i < tab->m; ++i)
				if (igual(tab->A[i][j], 1))
					x[k++] = tab->A[i][0];
		}
	}
}

Tableau *auxiliar(Tableau *tab) {
	int i, j;
	Tableau *aux = (Tableau*) malloc(sizeof(Tableau));
	aux->m = tab->m;
	aux->n = tab->n;

	for (i = 0; i < tab->m; ++i)
		for (j = 0; j < tab->n; ++j)
			aux->A[i][j] = tab->A[i][j];

	for (i = 1; i < aux->m; ++i)
		if (aux->A[i][0] < 0)
			for (j = 0; j < aux->n; ++j)
				aux->A[i][j] = -aux->A[i][j];

	for (j = 0; j < tab->n; ++j) {
		double soma = 0;
		for (i = 1; i < aux->m; ++i) {
			soma -= aux->A[i][j];
		}
		aux->A[0][j] = soma;
	}

	return aux;
}

int tem_matriz_identidade(Tableau *tab) {
	int i, j;
	for (i = 1; i < tab->m; i++) {
		for (j = tab->n - tab->m + 1; j < tab->n; j++) {
			printf("A[%d][%d] ", i + tab->n - tab->m, j);
//			printf("%lf ",tab->A[i][j]);
//			if(!igual(tab->A[i][j], 0)) printf("%lf ",tab->A[i][j]);
//			if (!igual(fabs(tab->A[i][j]), 1) && !igual(tab->A[j][i+tab->n - tab->m], 0)){
//				printf("[%lf %lf]",tab->A[i][j], tab->A[j][i+tab->n - tab->m]);
//				return 0;
//			}
			if (i + tab->n - tab->m == j && !igual(fabs(tab->A[i][j]), 1))
				return 0;
			else if (i + tab->n - tab->m != j && !igual(tab->A[i][j], 0))
				return 0;
		}
		printf("\n");
	}
	return 1;
}

int simplex(int m, int n, double z[], double A[][MAX], double b[], double *z0,
		double x[]) {
	int i, j, lin_pivo, col_pivo;
//Atenção o simplex foi invocado e a entrada vai ser mostrada
	/*	 printf("m = %d, n = %d \n",m,n);
	 mostra_A(A, m, n);
	 mostra_z(z, n);
	 mostra_b(b, m);*/

//TODO verificar se todas as colunas são linearmente independentes?
//Controi o tableau
	Tableau *tab = cria_tableau(m, n, z, A, b);
	mostra_tableau(tab);

	if (tem_matriz_identidade(tab))
		printf("\ntem matriz identidade\n");
	else
		printf("\nNÃOOOO TEM identidade\n");

	if (!tem_matriz_identidade(tab)) {
		// PROGRAMA AUXILIAR
		Tableau *aux = auxiliar(tab);
		mostra_tableau(aux);

		add_folgas(aux);
		printf("\nFolgas adicionadas\n");
		mostra_tableau(aux);
		//simplex
		while (coluna_pivo(aux) > -1) {
			//Acha a coluna pivô
			col_pivo = coluna_pivo(aux);
			lin_pivo = linha_pivo(aux, col_pivo); //Acha a linha pivô
			if (lin_pivo >= 0)
				pivoteamento(aux, lin_pivo, col_pivo); //Computa o novo tableau
			else
				return UNBD;
		}

		//Todos os valores na linha 0 são não-negativos
		if (igual(aux->A[0][0], 0)) {
			for (i = 1; i < tab->m; ++i) {
				for (j = 0; j < tab->n; ++j) {
					tab->A[i][j] = aux->A[i][j];
				}
			}
			mostra_tableau(tab);
		} else if (aux->A[0][0] > 0) {
			return UNBD;
		} else {
			return NFEA;
		}

	}

// A base é viável?
	if (!base_viavel(b, m)) { // Não! vamos criar um programa auxiliar e resolvê-lo
		printf("\nVamos precisar de um programa auxiliar!\n");
		Tableau *aux = auxiliar(tab);
		mostra_tableau(aux);

		//simplex
		while (coluna_pivo(aux) > -1) {
			//Acha a coluna pivô
			col_pivo = coluna_pivo(aux);
			lin_pivo = linha_pivo(aux, col_pivo); //Acha a linha pivô
			if (lin_pivo >= 0)
				pivoteamento(aux, lin_pivo, col_pivo); //Computa o novo tableau
			else
				return UNBD;
		}

		//Todos os valores na linha 0 são não-negativos
		if (igual(aux->A[0][0], 0)) {
			for (i = 1; i < tab->m; ++i) {
				for (j = 0; j < tab->n; ++j) {
					tab->A[i][j] = aux->A[i][j];
				}
			}
			mostra_tableau(tab);
		} else if (aux->A[0][0] > 0) {
			return UNBD;
		} else {
			return NFEA;
		}

	}	// Deu boa, podemos fazer o pivoteamento
	printf("\nA base é viável, mãos a obra!\n");
	while (coluna_pivo(tab) > -1) {
		//Acha a coluna pivô
		col_pivo = coluna_pivo(tab);
		lin_pivo = linha_pivo(tab, col_pivo); //Acha a linha pivô
		if (lin_pivo >= 0)
			pivoteamento(tab, lin_pivo, col_pivo); //Computa o novo tableau
		else
			return UNBD;
	}
//Todos os valores na linha 0 são não-negativos
	if (m > n) {
		int cont = 0;
		for (i = 1; i < tab->m; ++i)
			if (!igual(tab->A[i][0], 0))
				cont++;
		if (cont > n)
			return NFEA;
	}
	*z0 = tab->A[0][0];
	vetor_otimo(tab, x);
	return FEA;
}
