#include <iostream>
#include "SLAE.h"
#include "Vector.h"
#include "BoundValProblem.h"

using namespace std;

int main()
{
   BoundValProblem bvp = BoundValProblem();

   bvp.form_elems("tests/test1/regions.txt");
   bvp.form_boundaries("tests/test1/boundary.txt"); 

   bvp.form_portrait();
   bvp.build_global_mat();

   bvp.first_bound();

   bvp.solve();

   bvp.print_results("tests/test1/results.txt");

}