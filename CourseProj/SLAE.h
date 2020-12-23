#pragma once
#include "Matrix.h"
using namespace std;

class SLAE
{
public:                    
   int maxiter;         // ������������ ���������� ��������
   real eps;            // �������� ��������� ������������� �������

   vector<real> pr;     // ������ ������ �����

   vector<real> t;      // ��������������� ������ ��� ���
   vector<real> tt;     // ��������������� ������ ��� ���
   vector<real> rk1;    // ������ ������� �� ����. �������� ���
   vector<real> zk1;    // ������ ������ �� ����. �������� ���
   vector<real> AtAzk1; // ��������������� ������ ��� ���

   SLAE(int N, int _maxiter, real _eps)
   {
      maxiter = _maxiter;
      eps = _eps;

      pr.resize(N);

      t.resize(N);
      tt.resize(N);
      rk1.resize(N);
      zk1.resize(N);
      AtAzk1.resize(N);
   }

   SLAE()
   {

   }

   // ����� ����������� ����������, ���������� ���������� ��������
   int conj_grad_method(vector<real>& xk1, vector<real>& res, Matrix& mat)
   {
      for (int i = 0; i < mat.N; i++)
         res[i] = xk1[i] = 0;

      mat.matrix_vector_mult(xk1, t, mat.ggl, mat.ggu);       // t = A * x0
      mat.matrix_vector_mult(pr - t, rk1, mat.ggu, mat.ggl);  // r0 = AT(f - A * x0)
      zk1 = rk1;

      int k = 1;
      while (k < maxiter)
      {
         mat.matrix_vector_mult(zk1, t, mat.ggl, mat.ggu);    // t = A * zk-1

         mat.matrix_vector_mult(t, AtAzk1, mat.ggu, mat.ggl); // AtAzk1 = At * A * zk-1
         real ak = (rk1 * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         real bk = rk1 * rk1;
         rk1 = rk1 - ak * AtAzk1;
         bk = (rk1 * rk1) / bk;
         zk1 = rk1 + bk * zk1;

         real disc = norm(rk1) / norm(pr); // ������������� �������

         if (disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }

   // ����� ����������� ���������� � ����������������� ��������
   // ������������� ��������, ���������� ���������� ��������
   int conj_grad_pred_method(vector<real>& xk1, vector<real>& res, Matrix& mat, SLAE& fac_slae, Matrix& fac_mat)
   {
      for(int i = 0; i < mat.N; i++)
         res[i] = xk1[i] = 0;

      mat.matrix_vector_mult(xk1, t, mat.ggl, mat.ggu);       // t = A * x0
      mat.matrix_vector_mult(pr - t, rk1, mat.ggu, mat.ggl);  // r0 = AT(f - A * x0)

      // ������ z0 = M-1 * r0
      fac_slae.pr = rk1;
      fac_slae.conj_grad_method(t, zk1, fac_mat);

      int k = 1;
      while(k < maxiter)
      {
         // ������ tt = M-1 * rk-1
         fac_slae.pr = rk1;
         fac_slae.conj_grad_method(t, tt, fac_mat);

         mat.matrix_vector_mult(zk1, t, mat.ggl, mat.ggu);    // t = A * zk-1
         mat.matrix_vector_mult(t, AtAzk1, mat.ggu, mat.ggl); // AtAzk1 = At * A * zk-1

         real ak = (tt * rk1) / (AtAzk1 * zk1);
         xk1 = xk1 + ak * zk1;
         real bk = tt * rk1;
         rk1 = rk1 - ak * AtAzk1;

         // ������ tt = M-1 * rk
         fac_slae.pr = rk1;
         fac_slae.conj_grad_method(t, tt, fac_mat);

         bk = (tt * rk1) / bk;
         zk1 = tt + bk * zk1;

         real disc = norm(rk1) / norm(pr);

         if(disc < eps)
            break;
         else
            k++;
      }

      res = xk1;
      return k;
   }

};

