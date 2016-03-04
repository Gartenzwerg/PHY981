#include<iostream>
#include<fstream>
#include<armadillo>
#include<string>
#include<cmath>

using namespace std;
using namespace arma;

ifstream fin;
double hw = 10.0; //h_bar*omega = 10 MeV
int nStates = 80; //total number of states
int N = 16; //total number of neutrons
int n[80], l[80];
float j[80], mj[80], tz[80];

void HartreeFock (mat h0, float ****v, vec &E, mat &C);

int main() {
	
	string buffer;
	int temp1, temp2, temp3, temp4, temp5, temp6;
	fin.open("spdata.dat");
	getline(fin, buffer);
	for(int i=0; i < nStates; i++) {
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
	
	mat h0 = zeros(nStates, nStates);
	for(int i=0; i < nStates; i++) {
		h0(i, i) = (2*n[i] + l[i] + 3.0/2.0)*hw;
	}


	vec E = zeros(nStates, 1);
	mat C = eye(nStates, nStates);
	
	HartreeFock(h0, v, E, C);

	return 0;
}

void HartreeFock (mat h0, float ****v, vec &E, mat &C) {

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
				for(int x = 0; x < 16; x++){
					densityMatrix(alpha, beta) += C(alpha, x)*C(beta, x);
				}
			}
		}
		pot = zeros(nStates, nStates);
		for(int alpha = 0; alpha < nStates; alpha++){
			for(int beta = alpha; beta < nStates; beta++){
				if(l[alpha] != l[beta] || j[alpha] != j[beta] || mj[alpha] != mj[beta]) {
					pot(alpha, beta) = 0;
					continue;
				}
				for(int gamma = 0; gamma < nStates; gamma++){
					for(int delta = 0; delta < nStates; delta++){
						pot(alpha, beta) += densityMatrix(gamma, delta)*v[alpha][gamma][beta][delta];
					}
				}
			pot(beta, alpha) = pot(alpha, beta);
			}	
		}
		h = h0 + pot;
		eig_sym(E, C, h);
		it++;
		diff = E - ePrev;
		if( abs(diff.max()) < threshold ) break;
		ePrev = E;
	}
	/*cout<<"************************************************************************\n";

	float energy = 0;
	for(int counter=0; counter<N; counter++) {
		cout<<counter<<"-th single-particle harmonic oscillator energy = "<<h0(counter, counter)<<"MeV."<<endl;
		energy += h0(counter, counter);
	}
	cout<<"	Total harmonic oscillator energy = "<<energy<<"MeV."<<endl;*/
	cout<<"************************************************************************\n";

	for(int counter=0; counter<40; counter++) {
		cout<<counter<<"-th single-particle Hatree-Fock energy = "<<E(counter)<<"MeV."<<endl;
	}

}
