/* rand example: guess the number */
#include <iostream>
#include <fstream>
#include "rand.h"      
#include "Matrix.h"
using namespace std;

double func(double rand_num){
  return 3.0 / 4.0 * (1 - rand_num * rand_num);
}

int main(){
  int N_rand = 1000000, n = 1000;
  double a = -1, b = 1;
  Matrix x(n), P(n);
  for (int i = 1; i <= n; i++){
    x(i) = (b - a) / (n - 1) * (i - 1) + a;
    P(i) = func(x(i));
  }

  double pmax = 1, xtry = 0, ptry = 0;
  Matrix A(N_rand);
  long seed1 = 1.0, seed2 = 100.0;
  bool test_switch = false;
  for (int i = 1; i <= N_rand; i++){
    test_switch = true;
    while (test_switch){ 
      double R1 = rand(seed1);
      xtry = a + (b - a) * R1;
      ptry = func(xtry);
      double R2 = rand(seed2);
      //cout << i << " " << R1 << " " << xtry << " " << ptry << " " << ptry / pmax << " " << R2 << endl;
      if (ptry / pmax >= R2){
        A(i) = xtry;
        test_switch = false;
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
