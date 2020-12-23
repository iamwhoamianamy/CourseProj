#pragma once
#include "SLAE.h"
#include "Test.h"
using namespace std;

class BoundValProblem
{
public:
   Matrix Global;               // ���������� �������
   Matrix Fac_global;           // �������� ������������ ���������� �������
   Matrix StiffMatrix;          // ��������� ������� ���������
   Matrix WeightMatrix;         // ��������� ������� ����
   Matrix G1, G2, M;            // ��������������� ������� ��� ���������� ���������
                                // ��������� ������ � �������� ������ ������

   SLAE Slae;                   // �������� ������� ��� ������������������
   SLAE Fac_slae;               // �������� ������� c �������������������

   vector<real> B;              // ���������� ������ ������ �����
   vector<real> LocB;           // ��������� ������ ������ �����
   vector<real> Solution;       // �������

   Test test;                   // ���������� � ��������� �������,
                                // ��������� ������ � ����� � ��������
                                // ��������

   int reg_count;                     // ���������� �����������
   int elem_count;                    // ���������� �������� ���������
   int node_count;                    // ���������� �����

   vector<real> Xw;                   // � ���������� ������ �����������
   vector<real> Yw;                   // Y ���������� ������ �����������
   vector<vector<int>> reg_board_i;   // ������� ������ �����������
   vector<vector<int>> bound_cond;    // ���������� � ��������� ��������
   vector<real> Xc;                   // ������������ ����� �� X
   vector<real> Yc;                   // ������������ ����� �� Y

   vector<real> Xn;                   // ������ ��������� � ��������� ��������
   vector<real> Yn;                   // ������ ��������� y ��������� ��������

   vector<real> f_in_elem;            // ������ �������� ������ ����� f
                                      // � ����� ��������� ��������

   vector<int> global_indices;        // ������ � �������� ����� � ���������� ����������
   vector<real> True;                 // ������ �������� ������� � �����

   real big_number = 2.4E+10;         // ������� ����� ��� ����� 1 �������� �������

   // ����������� ������
   BoundValProblem()
   {
      reg_count = 3;

      global_indices.resize(9);

      Xn.resize(3);
      Yn.resize(3);

      f_in_elem.resize(9);
      LocB.resize(9);
   }

   // ������ ����� � ������������� ������
   void form_elems(string file)
   {
      // ������ ��������������� ������ ��� ���������� ������ ��������� � �����
      G1 = Matrix(9);
      G1.read_di_ggl("G1.txt");

      G2 = Matrix(9);
      G2.read_di_ggl("G2.txt");

      M = Matrix(9);
      M.read_di_ggl("M1.txt");

      // ������ ���������� �� ��������
      ifstream fin;
      fin.open(file);

      int x_count;                       // ���������� ������ �� ��� X
      fin >> x_count;
      Xw.resize(x_count);

      // ������ ������ �� ��� X
      for (int i = 0; i < reg_count; i++)
         fin >> Xw[i];

      int y_count;                       // ���������� ������ �� ��� Y
      fin >> y_count;
      Yw.resize(y_count);

      // ������ ������ �� ��� Y
      for (int i = 0; i < y_count; i++)
         fin >> Yw[i];

      reg_board_i.resize(reg_count);

      // ���������� �������� ������ �����������
      for (int i = 0; i < reg_count; i++)
      {
         reg_board_i[i].resize(4);
         for (int j = 0; j < 4; j++)
            fin >> reg_board_i[i][j];
      }

      vector<real> t(reg_count - 1);
      int n_div = 1;

      // ���������� � ������������ ������������ ���� �� X
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

      // ���������� � ������������ ������������ ���� �� Y
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

      // ������� ���������� �������� ���������
      elem_count = (Xc.size() - 1) * (Yc.size() - 1);

      // ������� ���������� ����� � ���������� ���������
      calc_global_indices(elem_count - 1);
      node_count = global_indices[8] + 1;

      
      True.resize(node_count);
      Slae = SLAE(node_count, 10000, 1e-14);
      Fac_slae = SLAE(node_count, 10000, 1e-14);

      Global.ig.resize(node_count + 1);
      B.resize(node_count);
      Solution.resize(node_count);

      Fac_global = Matrix(node_count, 0);
   }

   // ������ ���������� � ������� ��������
   void form_boundaries(string file)
   {
      ifstream fin;
      fin.open(file);

      int n;
      fin >> n;

      bound_cond.resize(n);

      for(int i = 0; i < n; i++)
      {
         bound_cond[i].resize(6);

         for(int j = 0; j < 6; j++)
            fin >> bound_cond[i][j];

      }
      fin.close();
   }

   // ��������� ��������� ������
   void gen_loc_matrices(real hx, real hy, real lam, real gam)
   {
      StiffMatrix = (lam / 90) * (hy / hx * G1 + hx / hy * G2);
      WeightMatrix = gam * hx * hy / 900 * M;
   }

   // ��������� ���������� ������� ������ �����
   void gen_loc_b(real hx, real hy, int reg_index)
   {
      for(int j = 0; j < 3; j++)
         for(int i = 0; i < 3; i++)
            f_in_elem[j * 3 + i] = test.f(Xn[i], Yn[j])[reg_index];

      M.matrix_vector_mult(f_in_elem, LocB, M.ggl, M.ggu);
      LocB = hx * hy / 900 * LocB;
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
   int get_reg_index()
   {
      real x = Xn[1];
      real y = Yn[1];
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
   void incert_to_row(int r, int c)
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
               incert_to_row(help[i][j][0], help[i][j][1]);
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

         int reg_index = get_reg_index();

         gen_loc_matrices(hx, hy, test.lambda()[reg_index], test.gamma()[reg_index]);

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

            for(int j = 0; j < 3; j++)
               for(int i = 0; i < 3; i++)
                  True[global_indices[j * 3 + i]] = test.u(Xn[i], Yn[j])[reg_index];
         }

         gen_loc_b(hx, hy, reg_index);

         for(int i = 0; i < 9; i++)
         {
            Global.di[global_indices[i]] += StiffMatrix.di[i] + WeightMatrix.di[i];
            B[global_indices[i]] += LocB[i];
         }

         int asdasd = 1;
      }
   }

   // ������� ����� ������ ������� �������
   void first_bound()
   {
      for(int elem_index = 0; elem_index < elem_count; elem_index++)
      {
         calc_node_coords(elem_index);
         calc_global_indices(elem_index);

         int reg_index = get_reg_index();

         for(int k = 0; k < bound_cond.size(); k++)
            if(bound_cond[k][0] == 1)
               for(int j = 0; j < 3; j++)
                  for(int i = 0; i < 3; i++)
                     if(Xn[i] >= Xw[bound_cond[k][2]] && Xn[i] <= Xw[bound_cond[k][3]] &&
                        Yn[j] >= Yw[bound_cond[k][4]] && Yn[j] <= Yw[bound_cond[k][5]])
                     {
                        Global.di[global_indices[j * 3 + i]] = big_number;
                        B[global_indices[j * 3 + i]] = big_number * test.ug(Xn[i], Yn[j])[bound_cond[k][1]];
                     }
      }
   }

   // ���������� �������
   void solve()
   {
      Slae.pr = B;
      Global.diag_fact(Fac_global);

      vector<real> x0(node_count, 0);
      cout << Slae.conj_grad_pred_method(x0, Solution, Global, Fac_slae, Fac_global) << endl;

      // ����� ����������� � �������
      for(int i = 0; i < node_count; i++)
      {
         cout << setw(3) << i + 1;
         cout << setw(12) << Solution[i];
         cout << setw(12) << True[i];
         cout << setw(15) << scientific << abs(Solution[i] - True[i]) << fixed << endl;
      }

      print_results("tests/test1/results.txt");
   }

   // ����� ����������� � ����
   void print_results(string file_name)
   {
      ofstream fout;
      fout.open(file_name);


      for(int i = 0; i < node_count; i++)
      {
         fout << i + 1;
         fout << "\t" << Solution[i];
         fout << "\t" << True[i];
         fout << "\t" << scientific << abs(Solution[i] - True[i]) << fixed << endl;
      }

      fout.close();
   }
};