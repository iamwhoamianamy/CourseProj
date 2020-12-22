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
   bvp.form_boundaries("tests/test1/boundary.txt"); 

   bvp.form_portrait();
   bvp.build_global_mat();
   bvp.first_bound();

   bvp.solve();

   cout << "Hello World!";
}