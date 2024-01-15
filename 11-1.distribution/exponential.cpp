/* rand example: guess the number */
#include <iostream>
#include <fstream>
#include "rand.h"      
#include "Matrix.h"
using namespace std;

double exponential_distribution(double rand_num){
  double lambda = 1.1;
  return exp( - rand_num / lambda) / lambda;
}

double func(double rand_num){
  return 3.0 / 4.0 * (1 - rand_num * rand_num);
}

int main(){
  int N_rand = 1000000, n = 100;
  double a = -1, b = 1;
  double h = (b - a) / (n - 1);
  Matrix A(N_rand), x(n), P(n), F(n);
  F.set(0.0);
  long seed = 10.0;
  for (int i = 1; i <= n; i++){
    x(i) = h * (i - 1) + a;
    //P(i) = exponential_distribution(x(i)); // [0, \intfy)
    P(i) = func(x(i)); // [-1, 1]
    for (int j = 1; j <= i; j++){
      F(i) += P(j) * h;
    }
   // cout << i << " " << x(i) << " " << P(i) << " " << F(i) << endl;
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
  ofstream randOut("rand.txt");
  for (int i = 1; i <= N_rand; i++){
    randOut << A(i) << endl;
  }
}
