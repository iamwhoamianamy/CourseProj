#pragma once
#include <vector>
using namespace std;
typedef double real;

class Test
{
public:
   vector<real> f(real x, real y)
   {
      return  { -4, -4, -4 };
   }

   vector<real> lambda()
   {
      return { 1, 1, 1 };
   }

   vector<real> gamma()
   {
      return { 0, 0, 0 };
   }

   vector<real> ug(real x, real y)
   {
      return { x * x + y * y, x * x + y * y, x * x + y * y };
   }

   vector<real> u(real x, real y)
   {
      return { x * x + y * y,  x * x + y * y, x * x + y * y };
   }
};

