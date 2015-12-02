/*
author: markus
  date: 2015-12-02
*/

#include <iostream>
#include <cmath>

using namespace std;

void f(double* y, double* k, double* k1, double* k2, double* a, const double eta, const double dx){

	k[0] = y[1] + dx*a[0]*k1[1] + dx*a[1]*k2[1];
	k[1] = (eta-pow(y[0] + dx*a[0]*k1[0] + dx*a[1]*k2[0],2)) * (y[0] + dx*a[0]*k1[0] + dx*a[1]*k2[0]);

}

int main(){

	const int N=10000;
	const double dx=100.0/N*1.0;
	const double eta=0.5;
	double y[2];
	double k1[2];
	double k2[2];
	double k3[2];
	double a[2];
	//double temp;

	y[0] = 1e-5;
	y[1] = sqrt(eta)*y[0];

	cout << 0 << "\t" << y[0] << "\t" << y[1] <<endl;

	for(int i=0; i<N; i++){

		a[0]=0; a[1]=0;
		f(y, k1, k1, k2, a, eta, dx);

		a[0]=0.5; a[1]=0;
		f(y, k2, k1, k2, a, eta, dx);

		a[0]=-1; a[1]=2;
		f(y, k3, k1, k2, a, eta, dx);

		y[0] = y[0] +dx/6*(k1[0]+4*k2[0]+k3[0]);
		y[1] = y[1] +dx/6*(k1[1]+4*k2[1]+k3[1]);

		cout << i*dx << "\t" << y[0] << "\t" << y[1] <<endl;
	}

	return 0;

}
