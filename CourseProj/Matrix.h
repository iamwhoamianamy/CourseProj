#pragma once
#include "Vector.h"
using namespace std;

class Matrix
{
public:
   int N;               // Размер матрицы
   int M;               // Количество элементов в треугольнике
   
   vector<int> ig;      // Указатели начала строк
   vector<int> jg;      // Номера столбцов внедиагональных элементов

   vector<real> ggl;    // Верхний треугольник
   vector<real> ggu;    // Нижний треугольник

   vector<real> di;     // Диагональ

   Matrix(string path)
   {
      ifstream fin;
      fin.open(path + "kuslau.txt");
      fin >> N;
      fin.close();

      read_vector<int>(path + "ig.txt", ig, N + 1);
      M = ig[N];
      read_vector<int>(path + "jg.txt", jg, M);
      read_vector<real>(path + "ggl.txt", ggl, M);
      read_vector<real>(path + "ggu.txt", ggu, M);
      read_vector<real>(path + "di.txt", di, N);
   }

   Matrix(int _N, int _M)
   {
      N = _N;
      M = _M;

      ggl.resize(M);
      ggu.resize(M);
      jg.resize(M);
      di.resize(N);
      ig.resize(N + 1);
   }
   Matrix()
   {

   }

   // Получение диагональной факторизации матрицы
   void diag_fact(Matrix& fact)
   {
      fact.di = di;

      for (int i = 0; i < N + 1; i++)
         fact.ig[i] = 0;
   }

   // Функция умножения матрицы на вектор vec, результат в res
   void matrix_vector_mult(const vector<real>& vec, vector<real>& res,
      const vector<real>& bot_tr, const vector<real>& top_tr)
   {
      for (int i = 0; i < N; i++)
         res[i] = 0;

      for (int i = 0; i < N; i++)
      {
         res[i] += vec[i] * di[i];

         int prof_len = ig[i + 1] - ig[i];
         for (int k = 0; k < prof_len; k++)
         {
            int i_in_gg = ig[i] + k;
            int j = jg[i_in_gg];
            res[i] += vec[j] * bot_tr[i_in_gg];
            res[j] += vec[i] * top_tr[i_in_gg];
         }
      }
   }

};