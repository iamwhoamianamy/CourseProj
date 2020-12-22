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

   Matrix(int _N)
   {
      N = _N;
      M = N * (N - 1) / 2;

      ggl.resize(M);
      ggu.resize(M);
      jg.resize(M);
      di.resize(N);
      ig.resize(N + 1);

      ig[0] = ig[1] = 0;

      di[0] = 1.0;

      for (int i = 1; i < N; i++)
      {
         int i0 = ig[i + 0];
         int i1 = ig[i + 1] = ig[i + 0] + i;

         for (int j = 0, k = i0; j < i; j++, k++)
            jg[k] = j;

         di[i] = 1 / (real(i) * 2 + 1);
      }
   }

   Matrix(const Matrix& mat)
   {
      N = mat.N;
      M = mat.M;

      ggl = mat.ggl;
      ggu = mat.ggu;
      di = mat.di;
      ig = mat.ig;
      jg = mat.jg;
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

   // Функция считывания диагонали и треугольников из файлов
   void read_di_ggl(string file_name)
   {
      ifstream fin;
      fin.open(file_name);

      for(int i = 0; i < N; i++)
         fin >> di[i];

      for(int i = 0; i < M; i++)
      {
         real t;
         fin >> t;
         ggl[i] = ggu[i] = t;
      }

      fin.close();
   }

   //// Функция вычисления суммы модулей всех элементов матрицы
   //real mat_sum_abs()
   //{
   //   real res = 0;
   //   res += vec_sum_abs(ggl);
   //   res += vec_sum_abs(ggu);
   //   res += vec_sum_abs(di);
   //   return res;
   //}
};

// Умножение матрицы на число
Matrix operator * (real val, const Matrix& mat)
{
   Matrix res = Matrix(mat);

   res.di = val * res.di;
   res.ggl = val * res.ggl;
   res.ggu = val * res.ggu;

   return res;
}

// Сложение матриц
Matrix operator + (const Matrix& mat1, const Matrix& mat2)
{
   Matrix res = Matrix(mat1);

   res.di = res.di + mat2.di;
   res.ggl = res.ggl + mat2.ggl;
   res.ggu = res.ggu + mat2.ggu;

   return res;
}