#pragma once
#include <vector>
using namespace std;
typedef double real;

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { 2 * (sin(x) + sin(y)),
//                2 * (sin(x) + sin(y)),
//                2 * (sin(x) + sin(y)) };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { sin(x) + sin(y),
//               sin(x) + sin(y),
//               sin(x) + sin(y) };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { sin(x) + sin(y),
//               sin(x) + sin(y),
//               sin(x) + sin(y) };
//   }
//};

class Test
{
public:
   vector<real> f(real x, real y)
   {
      return  { -4 + 1 * (x * x + y * y),
                -4 + 1 * (x * x + y * y),
                -4 + 1 * (x * x + y * y) };
   }

   vector<real> lambda()
   {
      return { 1,
               1,
               1 };
   }

   vector<real> gamma()
   {
      return { 1,
               1,
               1 };
   }

   vector<real> ug(real x, real y)
   {
      return { x * x + y * y,
               x * x + y * y,
               x * x + y * y };
   }

   vector<real> u(real x, real y)
   {
      return { x * x + y * y,
               x * x + y * y,
               x * x + y * y };
   }
};

//class Test
//{
//public:
//   vector<real> f(real x, real y)
//   {
//      return  { x + y,
//                x + y,
//                x + y };
//   }
//
//   vector<real> lambda()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> gamma()
//   {
//      return { 1, 1, 1 };
//   }
//
//   vector<real> ug(real x, real y)
//   {
//      return { x + y,
//               x + y,
//               x + y };
//   }
//
//   vector<real> u(real x, real y)
//   {
//      return { x + y,
//               x + y,
//               x + y };
//   }
//};
