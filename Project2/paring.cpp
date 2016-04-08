#include <iostream>
#include <armadillo>
#include <string>
#include <cmath>

using namespace std;
using namespace arma;

//define Slater determinants consisting of no more than two pairs of nucleons
struct SlaterDeterminant {
	// 'first' and 'second' store the p-value(s) of the pair(s); -1 indicates null.
	int first = -1;
	int second = -1;
};

struct TwoBodyOperator {
	int creation = -1;
	int annihilation = -1;
};

void Setup(SlaterDeterminant SD[], int N, TwoBodyOperator TB[], int NP, int pMax);

int main(int argc, char* argv[]) { 
	
	int N = 0, pMax = 0, NP; // N: number of Slater determinants; pMax: max p of occupied states; NP: number of TwoBodyOperator
	if(argv[1][0] == 'b') {
		N = 2;
		pMax = 2;
		NP = 4;
	} else if(argv[1][0] == 'c') {
		N = 6;
		pMax = 4;
		NP = 16;
	}

	double g = atof(argv[2]);

	SlaterDeterminant *SD;
	SD = new SlaterDeterminant [N];
	TwoBodyOperator *TB;
	TB = new TwoBodyOperator [NP];
	//set up the Slater determinants and two-body operators
	Setup(SD, N, TB, NP, pMax);

	/*for(int i = 0; i < N; i++) {
		cout<<SD[i].first<<"	"<<SD[i].second<<endl;
	}*/

	mat h = zeros(N, N);
	//set up the diagonal elements
	for(int i = 0; i < N; i++) {
		h(i, i) += 2*(SD[i].first - 1);
		if(SD[i].second != -1) h(i, i) += 2*(SD[i].second - 1);
		//cout<<h(i, i)<<endl;
	}

	for(int i = 0; i < N; i++) {
		for(int j = i; j < N; j++) {
			int temp = 0;
			SlaterDeterminant copy;
			for(int k = 0; k < NP; k++) {
				copy = SD[j];
				if(TB[k].annihilation == SD[j].first) {
					if((TB[k].creation < SD[j].second) || (SD[j].second == -1)) {
						copy.first = TB[k].creation;
					} else if(TB[k].creation > SD[j].second) {
						copy.first = SD[j].second;
						copy.second = TB[k].creation;
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second)) temp++;
					continue;
				}
				if((TB[k].annihilation == SD[j].second) && (SD[j].second != -1)) {
					if(TB[k].creation > SD[j].first) {
						copy.second = TB[k].creation;
					} else if(TB[k].creation < SD[j].first) {
						copy.second = SD[j].first;
						copy.first = TB[k].creation;
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second)) temp++;
				}
			}
			h(i, j) -= g*temp;
			h(j, i) = h(i, j);
			//cout<<h(i, j)<<"	";
		}
		//cout<<endl;
	}

	vec E = eig_sym(h);
	for(int i = 0; i < N; i++) {
		cout<<E(i)<<endl;
	}

	return 0;
}

void Setup(SlaterDeterminant SD[], int N, TwoBodyOperator TB[], int NP, int pMax) {
	if(N == 2) {
		SD[0].first = 1;
		SD[1].first = 2;
	} else if(N == 6) {
		int counter = 0;
		int upperLimit;
		for(int j = 1; j < 4; j++){
			upperLimit = counter + (4 - j);
			for(int i = counter; i < upperLimit; i++){
				SD[i].first = j;
				SD[i].second = (i-counter)+j+1;
			}
			counter += (4-j);
		}
	}
	int counter = 0;
	for(int i = 1; i < pMax+1; i++) {
		for(int j = 1; j < pMax+1; j++) {
			TB[counter].creation = i;
			TB[counter].annihilation = j;
			counter ++;
		}
	}
}