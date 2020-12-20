#include<bits/stdc++.h>
#include"simplex.h"

static const double eps = 1e-9;

typedef struct tableau {
	int m, n;
	double A[MAX][MAX];
} Tableau;

int igual(double x, double y) {
	return fabs(x - y) < eps;
}

int b_naonegativo(double b[], int m) {
	int i;
	for (i = 1; i <= m; ++i)
		if (b[i] < 0)
			return 0;
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

// Encontra a coluna pivô (variável que entra na base)
int coluna_pivo(Tableau *tab) {
	int i;
	for (i = 1; i < tab->n; ++i)
		if (tab->A[0][i] + eps < 0)
			return i;
	return -1;
}

//Encontra a linha pivô (variável que sai da base)
int linha_pivo(Tableau *tab, int col_pivo) {
	double razao;
	int i, lin_pivo;
	double min = INFINITY;

	for (i = 1; i < tab->m; ++i) {
		if (tab->A[i][col_pivo] - eps > 0) {
			razao = tab->A[i][0] / tab->A[i][col_pivo];
			if (min > razao) {
				min = razao;
				lin_pivo = i;
			}
		}
	}
	//Se não tem elemento positivo na coluna pivô
	if (min == INFINITY)
		return -1;

	return lin_pivo;
}

//Computa o novo tableau dado o pivô A[lin_pivo][lin_pivo]
void novo_tableau(Tableau *tab, int lin_pivo, int col_pivo) {
	int i, j;
	double pivo;
	pivo = tab->A[lin_pivo][col_pivo];

	//Nova linha pivô
	for (j = 0; j < tab->n; ++j) {
		tab->A[lin_pivo][j] = tab->A[lin_pivo][j] / pivo;
	}
	// Calcula as novas linhas
	for (i = 0; i < tab->m; ++i) {
		double tmp = tab->A[i][col_pivo];
		if (i != lin_pivo) {
			for (j = 0; j < tab->n; ++j) {
				tab->A[i][j] -= tab->A[lin_pivo][j] * tmp;
			}
		}
	}
}

/*Usa a regra de Bland para definir as variáveis que entram e que saem da base.
 Devolve FEA se programa é viável e UNBD se for ilimitado*/
int pivoteamento(Tableau *tab) {
	int col_pivo, lin_pivo;
	//Dado que com a regra de bland o simplex não cicla, o laço pára
	while (1) {
		col_pivo = coluna_pivo(tab); //Acha a coluna pivô
		if (col_pivo == -1)
			break;
		lin_pivo = linha_pivo(tab, col_pivo); //Acha a linha pivô
		if (lin_pivo > -1)
			novo_tableau(tab, lin_pivo, col_pivo); //Computa o novo tableau
		else
			return UNBD;
	}
	return FEA;
}

//Salva os valores encontrados para o vetor x[]
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
/* Devolve o tableau de um programa auxiliar. Caso existam valores negativos em b[],
 são invertidos os sinais dos elementos da restrição i cujo b_i < 0 */
Tableau *auxiliar(Tableau *tab) {
	int i, j;
	Tableau *aux = (Tableau*) malloc(sizeof(Tableau));
	aux->m = tab->m;
	aux->n = tab->n;

	for (i = 0; i < tab->m; ++i) {
		for (j = 0; j < tab->n; ++j) {
			//Inverte os sinais dos elementos da restrição i cujo b_i seja negativo
			if (tab->A[i][0] + eps < 0)
				aux->A[i][j] = -tab->A[i][j];
			else
				aux->A[i][j] = tab->A[i][j];
		}
	}

	for (j = 0; j < tab->n; ++j) {
		double soma = 0;
		for (i = 1; i < aux->m; ++i)
			soma -= aux->A[i][j];
		aux->A[0][j] = soma;
	}
	return aux;
}

/*Verifica se o programa possui base viável, devolve 1 se a submatriz das ultimas m
 colunas de A é uma matriz identidade e zero caso contrario*/
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

/* Verifica se existe variavel artificial na base. Se algum elemento de b[] for zero, procura por alguma variavel não nula do programa original, senão encontra,
mantém a variavel auxiliar igual a zero.*/
void var_artificial_na_base(Tableau *aux) {
	int i, j;
	for (i = 1; i < aux->m; ++i) {
		if (igual(aux->A[i][0], 0)) {
			for (j = 1; j < aux->n - aux->m + 1; ++j)
				if (!igual(aux->A[i][j], 0))
					novo_tableau(aux, i, j);
		}
	}
}

int simplex(int m, int n, double z[], double A[][MAX], double b[], double *z0,
		double x[]) {
	int i, j, lin_pivo, col_pivo;

	//Constroi o tableau
	Tableau *tab = cria_tableau(m, n, z, A, b);

	//Se a base não é viável, cria variaveis artificiais e um programa auxiliar
	if (!tem_matriz_identidade(tab)) {
		// PROGRAMA AUXILIAR
		Tableau *aux = auxiliar(tab);
		add_folgas(aux);

		//Faz pivoteamento enquanto existir variavel não negativa na função objetivo
		if (pivoteamento(aux))
			return UNBD;

		//Otimo do programa auxiliar encontrado
		if (igual(aux->A[0][0], 0)) {
			//se tem variavel artificial na base faz pivoteamento quando necessario
			var_artificial_na_base(aux);
			for (i = 1; i < tab->m; ++i)
				for (j = 0; j < tab->n; ++j)
					tab->A[i][j] = aux->A[i][j];
			free(aux);
		} else if (aux->A[0][0] - eps > 0) {
			free(aux);
			free(tab);
			return UNBD;
		} else {
			free(aux);
			free(tab);
			return NFEA;
		}
	}

	// A base é viável?
	if (!b_naonegativo(b, m)) { // Não! cria um programa auxiliar e resolve-o
		Tableau *aux = auxiliar(tab);

		if (pivoteamento(aux))
			return UNBD;

		if (igual(aux->A[0][0], 0)) {
			for (i = 1; i < tab->m; ++i)
				for (j = 0; j < tab->n; ++j)
					tab->A[i][j] = aux->A[i][j];
			free(aux);
		} else if (aux->A[0][0] - eps > 0) {
			free(aux);
			free(tab);
			return UNBD;
		} else {
			free(aux);
			free(tab);
			return NFEA;
		}
	}

	if (pivoteamento(tab)) {
		free(tab);
		return UNBD;
	} else { // otimo encontrado
		*z0 = tab->A[0][0];
		vetor_otimo(tab, x);
		free(tab);
		return FEA;
	}
}
