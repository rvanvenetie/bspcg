#include <cstdlib> 
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

typedef struct mat{
	int n;
	vector<int> i,j;
	vector<double> val;
} mat;

//Random number [0,1]
double rand_double(void) {
  return ( (double) rand() / (double) RAND_MAX);
}

mat rand_mat( int n,double sparsity) {
	vector<double> diag(n, 0.0);
	int cnt = 0;
	mat result;
	result.n = n;
	for (int i =1; i <= n; i++)
		for (int j=1; j < i; j++)
			if (rand_double() < sparsity) {
			  result.i.push_back(i);
				result.j.push_back(j);
				double new_val = 2*rand_double() -1;
				result.val.push_back(new_val);
				diag[i] += abs(new_val);
				diag[j] += abs(new_val);
			}
	double max = *max_element(diag.begin(), diag.end());
	//Loop over diagonal

	for (int i = 1; i <= n; i++)
	{
		result.i.push_back(i);
		result.j.push_back(i);
		double new_val = max;
		if (rand_double() < sparsity) 
			new_val += rand_double() * 2;
		result.val.push_back(new_val);
	}
	return result;
}

void print_mat(mat print) {
	printf("%%%%MatrixMarket matrix coordinate real symmetric\n");
	printf("%d %d %d\n", print.n, print.n, print.i.size());
	for (int i = 0; i < print.i.size(); i++)
		printf("%d %d %g\n", print.i[i], print.j[i], print.val[i]);
}
int main(int argc, char *argv[]) {
	if (argc != 3)
	{
		printf("Parameters are <n> <density>\n");
		printf("For example 50 0.01\n");
		return 0;
	}
	int n = atoi(argv[1]);
	double sparsity = atof(argv[2]);

	//Give every combination n and sparsity a different seed, so we have constant matrix per combination
	srand(n);
	srand((int) ((double) rand() / sparsity));
  print_mat(rand_mat(n, sparsity));

	return 0;
}
