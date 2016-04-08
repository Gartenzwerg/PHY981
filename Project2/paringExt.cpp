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
	int third = -1;
	int fourth = -1;
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
	} else if(argv[1][0] == '6') {
		N = 20;
		pMax = 6;
		NP = 36;
	} else if(argv[1][0] == '8') {
		N = 70;
		pMax = 8;
		NP = 64;
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
		if(SD[i].third != -1) h(i, i) += 2*(SD[i].third - 1);
		if(SD[i].fourth != -1) h(i, i) += 2*(SD[i].fourth - 1);
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
					} else if( TB[k].creation > SD[j].second ) {
						copy.first = SD[j].second;
						copy.second = TB[k].creation;
						if((SD[j].third != -1) && (TB[k].creation > SD[j].third)) {
							copy.second = SD[j].third;
							copy.third = TB[k].creation;
							if((SD[j].fourth != -1) && (TB[k].creation > SD[j].fourth)) {
								copy.third = SD[j].fourth;
								copy.fourth = TB[k].creation;
							}
						}
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second) && (SD[i].third == copy.third) && (SD[i].fourth == copy.fourth)) temp++;
					continue;
				}
				if((TB[k].annihilation == SD[j].second) && (SD[j].second != -1)) {
					if(TB[k].creation > SD[j].first) {
						copy.second = TB[k].creation;
						if((SD[j].third != -1) && (TB[k].creation > SD[j].third)) {
							copy.second = SD[j].third;
							copy.third = TB[k].creation;
							if((SD[j].fourth != -1) && (TB[k].creation > SD[j].fourth)) {
								copy.third = SD[j].fourth;
								copy.fourth = TB[k].creation;
							}
						}
					} else if(TB[k].creation < SD[j].first) {
						copy.second = SD[j].first;
						copy.first = TB[k].creation;
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second) && (SD[i].third == copy.third) && (SD[i].fourth == copy.fourth)) temp++;
					continue;
				}
				if((TB[k].annihilation == SD[j].third) && (SD[j].third != -1)) {
					if(TB[k].creation < SD[j].fourth || SD[j].fourth == -1) {
						copy.third = TB[k].creation;
						if(TB[k].creation < SD[j].second) {
							copy.second = TB[k].creation;
							copy.third = SD[j].second;
							if(TB[k].creation < SD[j].first) {
								copy.first = TB[k].creation;
								copy.second = SD[j].first;
							}
						}
					} else if((TB[k].creation > SD[j].fourth) && (SD[j].fourth != -1)) {
						copy.third = copy.fourth;
						copy.fourth = TB[k].creation;
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second) && (SD[i].third == copy.third) && (SD[i].fourth == copy.fourth)) temp++;
					continue;
				}
				if((TB[k].annihilation == SD[j].fourth) && (SD[j].fourth != -1)) {
					if(TB[k].creation > SD[j].third) {
						copy.fourth = TB[k].creation;
					} else if(TB[k].creation < SD[j].third) {
						copy.third = TB[k].creation;
						copy.fourth = SD[j].third;
						if(TB[k].creation < SD[j].second) {
							copy.second = TB[k].creation;
							copy.third = SD[j].second;
							if(TB[k].creation < SD[j].first) {
								copy.first = TB[k].creation;
								copy.second = SD[j].first;
							}
						}
					} else {
						continue;
					}
					if((SD[i].first == copy.first) && (SD[i].second == copy.second) && (SD[i].third == copy.third) && (SD[i].fourth == copy.fourth)) temp++;
					continue;
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
	} else if(N == 20) {
		int counter = 0;
		int upperLimit;
		for(int j = 1; j < pMax-1; j++) {
			for(int i = j + 1; i < pMax; i++) {
				upperLimit = counter + (pMax - i);
				for(int k = counter; k < upperLimit; k++) {
					SD[k].first = j;
					SD[k].second = i;
					SD[k].third = (k-counter)+i+1;
				}
				counter += (pMax-i);
			}
		}
	} else if(N == 70) {
		int counter = 0;
		int upperLimit;
		for(int j = 1; j < pMax-2; j++) {
			for(int i = j+1; i < pMax-1; i++) {
				for(int k = i+1; k < pMax; k++) {
					upperLimit = counter+(pMax - k);
					for(int l = counter; l < upperLimit; l++) {
						SD[l].first = j;
						SD[l].second = i;
						SD[l].third = k;
						SD[l].fourth = (l-counter)+k+1;
					}
					counter += (pMax - k);
				}
			}
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