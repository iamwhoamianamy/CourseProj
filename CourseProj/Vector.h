#pragma once
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

typedef double real;

// Умножение вектора на число
vector<real> operator * (real val, const vector<real>& vec)
{
   vector<real> res(vec.size());

   for (size_t i = 0; i < vec.size(); ++i)
      res[i] = val * vec[i];
   return res;
}

// Сложение векторов
vector<real> operator + (const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
      throw("a.size() != b.size()");

   vector<real> res(vec1.size());

   for (size_t i = 0; i < vec1.size(); ++i)
      res[i] = vec1[i] + vec2[i];
   return res;
}

// Вычитание векторов
vector<real> operator - (const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
     throw("a.size() != b.size()");

   vector<real> res(vec1.size());

   for (size_t i = 0; i < vec1.size(); ++i)
      res[i] = vec1[i] - vec2[i];
   return res;
}

// Скалярное произведение векторов
real operator *(const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
      throw("vec1.size() != vec2.size()");

   int n = vec1.size();
   real res = 0;

   for (int i = 0; i < n; i++)
      res += vec1[i] * vec2[i];

   return res;
}

// Норма вектора
real norm(const vector<real>& vec)
{
   return sqrt(vec * vec);
}

//// Функция вычисления суммы модулей всех элементов вектора
//real vec_sum_abs(vector<real>& vec)
//{
//   real res = 0;
//   for(int i = 0; i < vec.size(); i++)
//      res += abs(vec[i]);
//   return res;
//}

// Ввод вектора vec размерности n из файла с именем file_name
template<class type>
void read_vector(string file_name, vector<type>& vec, int n)
{
   vec.resize(n);
   ifstream fin;
   fin.open(file_name);
   for (int i = 0; i < n; i++)
      fin >> vec[i];

   fin.close();
}

// Вывод вектора vec в файл с именем file_name
template<class type>
void print_vector(string file_name, vector<type>& vec)
{
   ofstream fout;
   fout.open(file_name);

   int n = vec.size();
   fout << setprecision(16);

   for (int i = 0; i < n; i++)
      fout << vec[i] << " ";

   fout.close();
}
