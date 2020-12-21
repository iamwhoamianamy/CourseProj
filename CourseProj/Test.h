#pragma once
#include <vector>
using namespace std;
typedef double real;

class Test
{
public:
   vector<real> f(real x, real y)
   {
      return { x + y +1, 1, 0 };
   }

   vector<real> lambda()
   {
      return { 1, 2, 0};
   }

   vector<real> gamma()
   {
      return { 1, 2, 1 };
   }
};