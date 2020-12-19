#include <iostream>
#include "SLAE.h"
#include "Vector.h"
#include "BoundValProblem.h"

using namespace std;

int main()
{
   BoundValProblem bvp = BoundValProblem();
   Matrix M = bvp.G1 + bvp.G2;

   bvp.form_subregions("tests/test1/regions.txt");

   cout << "Hello World!";
}