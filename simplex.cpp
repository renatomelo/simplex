#include"simplex.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

static const double eps = 1.0e-8;

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

int coluna_pico(Tableau *tab) {
	int i;
	for (i = 1; i < tab->n; ++i)
		if (tab->A[0][i] < 0)
			return i;

	return -1;
}

//Encontra a linha pivô: menor razão p_i/x_v_i.
int linha_pivo(Tableau *tab, int col_pivo) {
	double razao;
	int i, lin_pivo;
	double min = INFINITY;

	for (i = 1; i < tab->m; ++i) {
		if (tab->A[i][col_pivo] > 0) {
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
		if (igual(tab->A[0][j], 0)) {
			for (i = 1; i < tab->m; ++i)
				if (igual(tab->A[i][j], 1))
					x[k++] = tab->A[i][0];
		}else
			x[k++] = 0;
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

	//Adiciona variaveis auxiliares em cada restrição
//	for (i = 1; i < aux->m; ++i) {
//		for (j = tab->n; j < aux->n; ++j) {
//			if ((j - tab->n) == i - 1)
//				aux->A[i][j] = 1;
//			else
//				aux->A[i][j] = 0;
//		}
//	}

	for (j = 0; j < tab->n; ++j) {
		double soma = 0;
		for (i = 1; i < aux->m; ++i) {
			soma -= aux->A[i][j];
		}
		aux->A[0][j] = soma;
	}

	return aux;
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

// Vai começar o simplex, prepare-se!
// Antes de tudo, a base é viável?
	if (!base_viavel(b, m)) { // Não! vamos criar um programa auxiliar e resolvê-lo
		printf("\nVamos precisar de um programa auxiliar!\n");
		Tableau *aux = auxiliar(tab);
		mostra_tableau(aux);

		//simplex
		for (j = 1; j < aux->n; ++j) {
			//Acha a coluna pivô
			if (aux->A[0][j] < 0) {
				col_pivo = j;
				lin_pivo = linha_pivo(aux, col_pivo); //Acha a linha pivô
				if (lin_pivo >= 0)
					pivoteamento(aux, lin_pivo, col_pivo); //Computa o novo tableau
				else
					return UNBD;
			}
		}

		if (j == aux->n) { //Todos os valores na linha 0 são não-negativos
			if (igual(aux->A[0][0], 0)) {
				for (j = 1; j < aux->n; ++j) {
					aux->A[0][j] = -z[j];
				}
				for (i = 0; i < tab->m; ++i) {
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

	}	// Deu boa, podemos fazer o pivoteamento
	printf("\nA base é viável, mãos a obra!\n");

	while (1) {
		//Acha a coluna pivô
		col_pivo = coluna_pico(tab);
		if (col_pivo > 0) {
			lin_pivo = linha_pivo(tab, col_pivo); //Acha a linha pivô
			if (lin_pivo >= 0)
				pivoteamento(tab, lin_pivo, col_pivo); //Computa o novo tableau
			else
				return UNBD;
		} else { //Todos os valores na linha 0 são não-negativos
			int cont = 0;
			if (m > n) {
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
	}
}
