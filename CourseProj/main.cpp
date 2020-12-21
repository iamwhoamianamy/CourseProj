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
   bvp.calc_global_indices(0);
   bvp.calc_node_coords(0);
   cout << bvp.get_reg_index(4, 2);

   bvp.form_portrait();
   bvp.build_global_mat();

   cout << "Hello World!";
}