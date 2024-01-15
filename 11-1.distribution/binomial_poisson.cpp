/* rand example: guess the number */
#include <iostream>
#include <fstream>
#include "rand.h"      
#include "Matrix.h"
using namespace std;

double factorial(int N){
  double tmp = 1;
  for (int i = 1; i <= N; i++){
    tmp *= i;
  }
  return tmp;
}

double binomial_distribution(int N_total, int N_count){
  return factorial(N_total) / factorial(N_count) / factorial(N_total - N_count) * pow(0.5, N_total);
}

double poisson_distribution(int N_count){
  double lambda = 10;
  return exp(-lambda) * pow(lambda, N_count) / factorial(N_count);
}

int main ()
{ 
  int N_rand = 1000000, n = 100;
  Matrix A(N_rand), x(n), P(n), F(n);
  F.set(0.0);
  long seed = 1.0;
  for (int i = 1; i <= n; i++){
    x(i) = i;
    P(i) = binomial_distribution(n, i);
    //P(i) = poisson_distribution(i);
    for (int j = 1; j <= i; j++){
      F(i) += P(j);
    }
//    cout << i << " " << x(i) << " " << P(i) << " " << F(i) << endl;
  }
  if (abs(F(n) - 1) > 1e-3){
    cout << "The cumulative distribution function is not complete!!" << endl;
    cout << "The summation of probability distribution function is " << F(n) << endl;
  }
  for (int i = 1; i <= N_rand; i++){
    double R = rand(seed);
    for (int j = 1; j <= n; j++){
      if (R < F(j)){
        A(i) = x(j);
        break;
      }
    }
  }
  ofstream randOut("rand.txt"), NOut("N.txt"), POut("P.txt");
  for (int i = 1; i <= n; i++){
    NOut << x(i) << endl;
    POut << P(i) << endl;
  }
  for (int i = 1; i <= N_rand; i++){
    randOut << A(i) << endl;
  }
}
