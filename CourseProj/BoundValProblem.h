#pragma once
#include "SLAE.h"
#include "Test.h"
using namespace std;

class BoundValProblem
{
public:
   Matrix Global;               // ���������� �������
   Matrix StiffMatrix;          // ��������� ������� ���������
   Matrix WeightMatrix;         // ��������� ������� ����
   Matrix G1, G2, M;            // ��������������� ������� ��� ���������� ���������
                                // ��������� ������ � �������� ������ ������
   vector<real> BLoc;           // ��������� ������ ������ �����

   Test test;                   // ���������� � ��������� �������,
                                // ��������� ������ � ����� � ��������
                                // ��������

   int reg_count;                     // ���������� �����������
   int elem_count;                    // ���������� �������� ���������
   int node_count;                    // ���������� �����

   vector<real> Xw;                   // � ���������� ������ �����������
   vector<real> Yw;                   // Y ���������� ������ �����������
   vector<vector<real>> reg_board_i;  // ������� ������ �����������
   vector<real> Xc;                   // ������������ ����� �� X
   vector<real> Yc;                   // ������������ ����� �� Y

   vector<real> Xn;                   // ������ ��������� � ��������� ��������
   vector<real> Yn;                   // ������ ��������� y ��������� ��������

   vector<real> lambda;
   vector<real> gamma;
   vector<real> F;

   vector<int> global_indices;        // ������ � �������� ����� � ���������� ����������

   // ������������� ������
   BoundValProblem()
   {
      G1 = Matrix(9);
      G1.read_di_ggl("G1.txt");

      G2 = Matrix(9);
      G2.read_di_ggl("G2.txt");

      M = Matrix(9);
      M.read_di_ggl("M1.txt");

      global_indices.resize(9);
      Xn.resize(3);
      Yn.resize(3);
   }

   // ������ �����
   void form_subregions(string file)
   {
      reg_count = 3;

      ifstream fin;
      fin.open(file);

      int x_count;
      fin >> x_count;
      Xw.resize(x_count);

      for (int i = 0; i < reg_count; i++)
         fin >> Xw[i];

      int y_count;
      fin >> y_count;
      Yw.resize(y_count);

      for (int i = 0; i < y_count; i++)
         fin >> Yw[i];

      reg_board_i.resize(reg_count);

      for (int i = 0; i < reg_count; i++)
      {
         reg_board_i[i].resize(4);
         for (int j = 0; j < 4; j++)
            fin >> reg_board_i[i][j];
      }

      vector<real> t(reg_count - 1);
      int n_div = 1;

      for (int i = 0; i < reg_count - 1; i++)
      {
         fin >> t[i];
         n_div += t[i];
      }

      fin.close();

      Xc.resize(n_div);

      int k = 0;
      for (int i = 0; i < reg_count - 1; i++)
      {
         Xc[k] = Xw[i];
         k++;
         for (int j = 0; j < t[i] - 1; j++, k++)
            Xc[k] = (Xw[i + 1] - Xw[i]) / t[i] + Xc[k - 1];
      }

      Xc[k] = Xw[x_count - 1];

      n_div = 1;

      for (int i = 0; i < reg_count - 1; i++)
      {
         fin >> t[i];
         n_div += t[i];
      }

      Yc.resize(n_div);

      k = 0;
      for (int i = 0; i < reg_count - 1; i++)
      {
         Yc[k] = Yw[i];
         k++;
         for (int j = 0; j < t[i] - 1; j++, k++)
         {
            Yc[k] = (Xw[i + 1] - Xw[i]) / t[i] + Yc[k - 1];
         }
      }
      Yc[k] = Yw[y_count - 1];

      elem_count = (Xc.size() - 1) * (Yc.size() - 1);
      calc_global_indices(elem_count - 1);
      node_count = global_indices[8] + 1;
   }

   // ��������� ��������� ������
   void gen_loc_matrices(real hx, real hy, real lam, real gam)
   {
      StiffMatrix = (lam / 90 ) * (hy / hx * G1 + hx / hy * G2);
      WeightMatrix = gam * hx * hy / 900 * M;
   }

   // ��������� ������ global_indices ���������, ���������������� ���������� ���������
   // ����� ��������� �������� � ������� elem_index(���������� � ����)
   void calc_global_indices(int elem_index)
   {
      elem_index++;
      int Nx = Xc.size();
      int k = 2 * floor((elem_index - 1) / (Nx - 1)) * (2 * Nx - 1) + 2 * ((elem_index - 1) % (Nx - 1) + 1);
      k -= 2;

      global_indices[0] = k + 0;
      global_indices[1] = k + 1;
      global_indices[2] = k + 2;

      global_indices[3] = k + 2 * Nx - 1;
      global_indices[4] = k + 2 * Nx;
      global_indices[5] = k + 2 * Nx + 1;

      global_indices[6] = k + 2 * (2 * Nx - 1);
      global_indices[7] = k + 2 * (2 * Nx - 1) + 1;
      global_indices[8] = k + 2 * (2 * Nx - 1) + 2;
   }

   // ������� ���������� ����� ��������� ��������
   // � ������� elem_index(���������� � ����)
   void calc_node_coords(int elem_index)
   {
      elem_index++;
      int Nx = Xc.size();

      int p = (elem_index - 1) % (Nx - 1) + 1;
      int s = floor((elem_index - 1) / (Nx - 1)) + 1;

      p--;
      s--;

      Xn[0] = Xc[p];
      Xn[1] = (Xc[p] + Xc[p + 1]) / 2;
      Xn[2] = Xc[p + 1];

      Yn[0] = Yc[s];
      Yn[1] = (Yc[s] + Yc[s + 1]) / 2;
      Yn[2] = Yc[s + 1];
   }

   // ������� ������ ���������� �� ����������� ������������ ���� ��������� ��������
   int get_reg_index(real x, real y)
   {
      for (int i = 0; i < reg_count; i++)
      {
         real left = Xw[reg_board_i[i][0]];
         real right = Xw[reg_board_i[i][1]];

         real bot = Yw[reg_board_i[i][2]];
         real top = Yw[reg_board_i[i][3]];

         if (x >= left && x <= right && y >= bot && y <= top)
            return i;
      }
   }

   // ��������������� ������� ��� ������������ ��������
   void incert_to_mat(int r, int c)
   {
      int i_in_jg = Global.ig[r];
      int prof_len = Global.ig[r + 1] - Global.ig[r];

      bool found = false;

      for(int k = i_in_jg; k < i_in_jg + prof_len; k++)
         if(Global.jg[k] == c)
         {
            found = true;
            break;
         }
      
      if(!found)
      {
         for(int l = r + 1; l < Global.ig.size(); l++)
            Global.ig[l]++;

         int k = i_in_jg;

         while((k < i_in_jg + prof_len) && Global.jg[k] < c)
            k++;

         Global.jg.insert(Global.jg.begin() + k, c);
      }
   }

   // ��������� ������� ��������� �������
   void form_portrait()
   {
      Global.ig.resize(node_count + 1);
      Global.jg;

      Global.ig[0] = Global.ig[1] = 0;

      for(int i = 0; i < elem_count; i++)
      {
         calc_global_indices(i);
         vector<vector<vector<int>>> help(9);

         for(int i = 0; i < 9; i++)
            help[i].resize(9);

         for(int i = 0; i < 9; i++)
            for(int j = 0; j < 9; j++)
               help[i][j] = { global_indices[i], global_indices[j] };


         for(int i = 1; i < 9; i++)
            for(int j = 0; j < i; j++)
               incert_to_mat(help[i][j][0], help[i][j][1]);
      }

      Global.N = Global.ig.size() - 1;
      Global.di.resize(Global.N);

      Global.M = Global.jg.size();
      Global.ggl.resize(Global.M);
      Global.ggu.resize(Global.M);
   }

   // ���������� �������� � �������
   void add_to_mat(int r, int c, real val_l, real val_u)
   {
      int beg_prof = Global.ig[r];
      int len_prof = Global.ig[r + 1] - Global.ig[r];

      for(int i_in_prof = beg_prof; i_in_prof < beg_prof + len_prof; i_in_prof++)
      {
         if(Global.jg[i_in_prof] == c)
         {
            Global.ggl[i_in_prof] += val_l;
            Global.ggu[i_in_prof] += val_u;
            break;
         }
      }
   }

   // ������ ���������� ������� 
   void build_global_mat()
   {
      for(int n_elem = 0; n_elem < elem_count; n_elem++)
      {
         calc_global_indices(n_elem);
         calc_node_coords(n_elem);

         real hx = Xn[2] - Xn[0];
         real hy = Yn[2] - Yn[0];

         int n_region = get_reg_index(Xn[1], Yn[1]);

         gen_loc_matrices(hx, hy, test.lambda()[n_region], test.gamma()[n_region]);

         for(int i = 1; i < 9; i++)
         {
            int beg_prof = StiffMatrix.ig[i];
            int len_prof = StiffMatrix.ig[i + 1] - StiffMatrix.ig[i];

            for(int i_in_prof = beg_prof; i_in_prof < beg_prof + len_prof; i_in_prof++)
            {
               int j = StiffMatrix.jg[i_in_prof];

               real val_l = StiffMatrix.ggl[i_in_prof] + WeightMatrix.ggl[i_in_prof];
               real val_u = StiffMatrix.ggu[i_in_prof] + WeightMatrix.ggu[i_in_prof];

               add_to_mat(global_indices[i], global_indices[j], val_l, val_u);
            }
         }

         int asdasd = 1;
      }


   }
};