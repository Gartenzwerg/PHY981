#include<iostream>
#include<fstream>
#include<armadillo>
#include<string>
#include<cmath>

using namespace std;
using namespace arma;

ifstream fin;
double hw = 10.0; //h_bar*omega = 10 MeV
int nStates = 40; //total number of states
int N = 8; //total number of neutrons
int n[80], l[80];
float j[80], mj[80], tz[80];

void HartreeFock (mat h0, float ****v, vec &E, mat &C);

int main() {
	
	string buffer;
	int temp1, temp2, temp3, temp4, temp5, temp6, counter=0;
	fin.open("spdata.dat");
	getline(fin, buffer);
	for(int i=0; i < 2*nStates; i++) {
		fin >> buffer >> buffer >> temp1 >> temp2 >> temp3 >> temp4 >> temp5 >> temp6;
		n[i] = temp2;
		l[i] = temp3;
		j[i] = temp4/2.0;
		mj[i] = temp5/2.0;
		tz[i] = temp6/2.0;
	}
	fin.close();

	float ****v = new float ***[80];
	for(int i=0; i<80; i++) {
		v[i] = new float **[80];
		for(int j=0; j<80; j++){
			v[i][j] = new float *[80];
			for(int k=0; k<80; k++){
				v[i][j][k] = new float [80];
			}
		}
	}
	
	fin.open("twobody.dat");
	getline(fin, buffer);
	getline(fin, buffer);
	int a, b, c, d;
	double temp;
	while(fin >> a >> b >> c >> d >> temp) {
		v[a-1][b-1][c-1][d-1] = temp;
	}
	fin.close();
	
	counter = 0;
	mat h0 = zeros(nStates, nStates);
	for(int i=0; i < 2*nStates; i++) {
		if(tz[i] > 0) {
			h0(counter, counter) = (2*n[i] + l[i] + 3.0/2.0)*hw;
			counter++;
		}
	}


	vec E = zeros(nStates, 1);
	mat C = eye(nStates, nStates);
	
	HartreeFock(h0, v, E, C);

	return 0;
}

void HartreeFock (mat h0, float ****v, vec &E, mat &C) {

	int counterA = 0, counterB = 0, counterC = 0, counterD = 0;
	mat densityMatrix = zeros(nStates, nStates), Hartree, pot, h;
	int iterationLimit = 100, it = 0;
	vec ePrev = zeros(nStates, 1);
	vec diff;
	double threshold = 1E-10;	

	while( it < iterationLimit ) {
		cout << "iteration = " << it <<endl;
		densityMatrix = zeros(nStates, nStates);
		for(int alpha = 0; alpha < nStates; alpha++){
			for(int beta = 0; beta < nStates; beta++){
				for(int x = 0; x < 8; x++){
					densityMatrix(alpha, beta) += C(alpha, x)*C(beta, x);
				}
			}
		}
		pot = zeros(nStates, nStates);
		for(int alpha = 0; alpha < 2*nStates; alpha++){
			if(tz[alpha] < 0) continue;
			counterB = counterA;
			for(int beta = alpha; beta < 2*nStates; beta++){
				if(tz[beta] < 0) continue;
				if(l[alpha] != l[beta] || j[alpha] != j[beta] || mj[alpha] != mj[beta]) {
					pot(counterA, counterB) = 0;
					counterB++;
					continue;
				}
				for(int gamma = 0; gamma < 2*nStates; gamma++){
					if(tz[gamma] < 0) continue;
					for(int delta = 0; delta < 2*nStates; delta++){
						if(tz[delta] < 0) continue;
						pot(counterA, counterB) += densityMatrix(counterC, counterD)*v[alpha][gamma][beta][delta];
						counterD++;
					}
				counterD = 0;
				counterC++;
				}
			pot(counterB, counterA) = pot(counterA, counterB);
			counterC = 0;
			counterB++;
			}
		counterA++;		
		}
		counterA = 0;
		h = h0 + pot;
		eig_sym(E, C, h);
		it++;
		diff = E - ePrev;
		if( abs(diff.max()) < threshold ) break;
		ePrev = E;
	}
	cout<<"************************************************************************\n";

	float energy = 0;
	for(int counter=0; counter<N; counter++) {
		cout<<counter<<"-th single-particle harmonic oscillator energy = "<<h0(counter, counter)<<"MeV."<<endl;
		energy += h0(counter, counter);
	}
	cout<<"	Total harmonic oscillator energy = "<<energy<<"MeV."<<endl;
	cout<<"************************************************************************\n";

	for(int counter=0; counter<N; counter++) {
		cout<<counter<<"-th single-particle Hatree-Fock energy = "<<E(counter)<<"MeV."<<endl;
	}

	energy = 0;
	for(int i=0; i < N; i++){
		for(int alpha = 0; alpha < nStates; alpha++){
			energy += C(alpha, i)*C(alpha, i)*h0(alpha, alpha);
		}
	}
	counterA = 0; counterB = 0; counterC = 0; counterD = 0;
	for(int m = 0; m < N; m++) {
		for(int n = 0; n < N; n++) {
			for(int alpha = 0; alpha < 2*nStates; alpha++){
				if(tz[alpha] < 0) continue;
				for(int beta = 0; beta < 2*nStates; beta++){
					if(tz[beta] < 0) continue;
					if(l[alpha] != l[beta] || j[alpha] != j[beta] || mj[alpha] != mj[beta]) {
						counterB++;
						continue;
					}
					for(int gamma = 0; gamma < 2*nStates; gamma++){
						if(tz[gamma] < 0) continue;
						for(int delta = 0; delta < 2*nStates; delta++){
							if(tz[delta] < 0) continue;
							energy += 0.5*C(counterA, m)*C(counterC, n)*C(counterB, m)*C(counterD, n)*v[alpha][gamma][beta][delta];
							counterD++;
						}
					counterD = 0;
					counterC++;
					}
				counterC = 0;
				counterB++;
				}
			counterB = 0;
			counterA++;		
			}
			counterA = 0;
			}
	}

	cout<<"	Total Hartree Fock Energy = "<<energy<<"MeV.\n";

}
