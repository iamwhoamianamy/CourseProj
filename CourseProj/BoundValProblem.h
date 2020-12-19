#pragma once
#include "SLAE.h"
using namespace std;

class BoundValProblem
{
public:
   Matrix Global;               // Глобальная матрица
   Matrix StiffMatrix;          // Локальная матрица жесткости
   Matrix WeightMatrix;         // Локальная матрица масс
   Matrix G1, G2, M;            // Вспомогательные матрицы для вычисления элементов
                                // локальных матриц и векторов правых частейы
   vector<real> BLoc;           // Локальный вектор правой части

   int subr_count;                  // Количество подобластей
   vector<real> Xw;                 // Х координаты границ подобластей
   vector<real> Yw;                 // Y координаты границ подобластей
   vector<vector<real>> regions_i;  // Индексы границ подобластей
   vector<real> coord_X;            // Координатные линии по X
   vector<real> coord_Y;            // Координатные линии по Y

   BoundValProblem()
   {
      G1 = Matrix(6);
      G1.di = { 28, 64, 28, 112, 256, 112 };
      G1.ggl = G1.ggu = {
         -32,
           4, -32,
          14, -16,   2,
         -16,  32, -16, -128,
           2, -16,  14,   16, -128 };

      G2 = Matrix(6);
      G2.di = { 28, 112, 28, 64, 256, 64 };
      G2.ggl = G2.ggu = {
          14,
          -7,   14,
         -32,  -16, 8,
         -16, -128, -16,  32,
           8,  -16, -32, -16, 32 };

      M = Matrix(6);
      M.di = { 16, 64, 16 ,64 ,256 , 64 };
      M.ggl = M.ggu = {
          8,
         -4,  8,
          8,  4, -2,
          4, 32,  4,  32,
         -2,  4,  8, -16, 32 };

   }

   // Генерация локальных матриц
   void gen_loc_matrices(real hx, real hy, real lam, real gam, real bet)
   {
      StiffMatrix = (lam / 90 ) * (hy / hx * G1 + hx / hy * G2);
      WeightMatrix = gam * hx * hy / 900 * M;
   }

   // Чтение сетки
   void form_subregions(string file)
   {
      subr_count = 3;

      ifstream fin;
      fin.open(file);

      int x_count;
      fin >> x_count;
      Xw.resize(x_count);

      for (int i = 0; i < subr_count; i++)
         fin >> Xw[i];

      int y_count;
      fin >> y_count;
      Yw.resize(y_count);

      for (int i = 0; i < y_count; i++)
         fin >> Yw[i];

      regions_i.resize(subr_count);

      for (int i = 0; i < subr_count; i++)
      {
         regions_i[i].resize(4);
         for (int j = 0; j < 4; j++)
            fin >> regions_i[i][j];
      }

      vector<real> t(subr_count - 1);
      int n_div = 1;

      for (int i = 0; i < subr_count - 1; i++)
      {
         fin >> t[i];
         n_div += t[i];
      }

      coord_X.resize(n_div);

      int k = 0;
      for (int i = 0; i < subr_count - 1; i++)
      {
         coord_X[k] = Xw[i];
         k++;
         for (int j = 0; j < t[i] - 1; j++, k++)
         {
            coord_X[k] = (Xw[i + 1] - Xw[i]) / t[i] + coord_X[k - 1];
         }
      }
      coord_X[k] = Xw[x_count - 1];

      n_div = 1;

      for (int i = 0; i < subr_count - 1; i++)
      {
         fin >> t[i];
         n_div += t[i];
      }

      coord_Y.resize(n_div);

      k = 0;
      for (int i = 0; i < subr_count - 1; i++)
      {
         coord_Y[k] = Yw[i];
         k++;
         for (int j = 0; j < t[i] - 1; j++, k++)
         {
            coord_Y[k] = (Xw[i + 1] - Xw[i]) / t[i] + coord_Y[k - 1];
         }
      }
      coord_Y[k] = Yw[y_count - 1];

      fin.close();
   }
};