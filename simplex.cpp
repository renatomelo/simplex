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

int b_positivo (double b[], int m) {
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

//Encontra a linha pivô
int linha_pivo(Tableau *tab, int col_pivo) {
	double razao;
	int i, lin_pivo;
	double min = INFINITY;

	for (i = 1; i < tab->m; ++i) {
		if (tab->A[i][col_pivo] - eps > 0) {
			razao = tab->A[i][0] / tab->A[i][col_pivo];
			if (min >= razao) {
				min = razao;
				lin_pivo = i;
			}
		}
	}
	//Se não tem elemento positivo na coluna pivô
	if (min == INFINITY) return -1;

	return lin_pivo;
}

//Usa regra de Bland para escolher a variavel que entra e a que sai
void pivoteamento(Tableau *tab, int lin_pivo, int col_pivo) {
	int i, j;
	double pivo;
	pivo = tab->A[lin_pivo][col_pivo];

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
		if (aux->A[i][0] + eps < 0)
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
			if (i + tab->n - tab->m == j && !igual(fabs(tab->A[i][j]), 1))
				return 0;
			else if (i + tab->n - tab->m != j && !igual(tab->A[i][j], 0))
				return 0;
		}
	}
	return 1;
}

int simplex(int m, int n, double z[], double A[][MAX], double b[], double *z0,
		double x[]) {
	int i, j, lin_pivo, col_pivo;

	//Controi o tableau
	Tableau *tab = cria_tableau(m, n, z, A, b);

	if (!tem_matriz_identidade(tab)) {
		// PROGRAMA AUXILIAR
		Tableau *aux = auxiliar(tab);
		add_folgas(aux);

		while (coluna_pivo(aux) > -1) {
			col_pivo = coluna_pivo(aux);
			lin_pivo = linha_pivo(aux, col_pivo);
			if (lin_pivo >= 0)
				pivoteamento(aux, lin_pivo, col_pivo); //Computa o novo tableau
			else
				return UNBD;
		}

		//Todos os valores na linha 0 são não-negativos
		if (igual(aux->A[0][0], 0)) {
			for (i = 1; i < tab->m; ++i)
				for (j = 0; j < tab->n; ++j)
					tab->A[i][j] = aux->A[i][j];
		} else if (aux->A[0][0] - eps > 0) {
			return UNBD;
		} else {
			return NFEA;
		}
	}

	// A base é viável?
	if (!b_positivo(b, m)) { // Não! vamos criar um programa auxiliar e resolvê-lo
		Tableau *aux = auxiliar(tab);

		while (coluna_pivo(aux) > -1) {
			col_pivo = coluna_pivo(aux);
			lin_pivo = linha_pivo(aux, col_pivo);
			if (lin_pivo >= 0)
				pivoteamento(aux, lin_pivo, col_pivo); //Computa o novo tableau
			else
				return UNBD;
		}

		//Todos os valores na linha 0 são não-negativos
		if (igual(aux->A[0][0], 0)) {
			for (i = 1; i < tab->m; ++i)
				for (j = 0; j < tab->n; ++j)
					tab->A[i][j] = aux->A[i][j];
		} else if (aux->A[0][0] - eps > 0) {
			return UNBD;
		} else {
			return NFEA;
		}
	}

	while (coluna_pivo(tab) > -1) {
		col_pivo = coluna_pivo(tab); //Acha a coluna pivô
		lin_pivo = linha_pivo(tab, col_pivo); //Acha a linha pivô
		if (lin_pivo > -1)
			pivoteamento(tab, lin_pivo, col_pivo); //Computa o novo tableau
		else
			return UNBD;
	}

	*z0 = tab->A[0][0];
	vetor_otimo(tab, x);
	return FEA;
}
