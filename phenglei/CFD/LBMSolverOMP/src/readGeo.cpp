#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include<sstream>
#include<algorithm>
#include "vector3.cpp"    //! ע��Ӧ����vector3.cpp���ǽ���vector3.hpp
using namespace std;

void getNumberFromString(string s, int *a, int n);
bool IntersectTriangle( vector3<double> &orig, vector3<double> &dir,
    vector3<double> &v0, vector3<double> &v1, vector3<double> &v2,
    double *t, double *u, double *v);

inline int index(int x, int y, int z, int NX, int NY) 
{
    return (z * NX * NY) + (y * NX) + x;
}


void read_Geo(int direction_size, vector3<int> *ei, int *obst, int *NN, int &num_node2, int *node2, double *delta)
{
    int NX = NN[0], NY = NN[1], NZ = NN[2];

    for (int k = 0; k < direction_size; k++)
    {
        cout << "(" << ei[k].x << " " << ei[k].y << "  " << ei[k].z << ")" << '\n';
    }

    char buffer[256];
    int ele_point = 0, total_ibpoint = 0, *a;

    double x_min, x_max, y_min, y_max, z_min, z_max;
    double x_length, y_length, z_length;
    double *lagx, *lagy, *lagz;
    int *Ep0, *Ep1, *Ep2;

    a = new int[2]();

    //! �����ļ�
    string filename1;
    filename1 = "mesh_3d.dat";
    if (direction_size == 9) filename1 = "mesh_2d.dat";
    ifstream infile(filename1);
    //ofstream outfile("out.plt");    //! ����ļ�
    if (!infile) 
    {
        cout << "Unable to open input complex geometry (shape) file:mesh_3d.dat or mesh_2d.dat";
        exit(1);    //! terminate with error
    
    }
    //if (!outfile) 
    //{
    //    cout << "Unable to open output file";
    //    exit(1);    // terminate with error
    //
    //}
    infile.getline(buffer, 60);    //! ����
    infile.getline(buffer, 60);    //! ����
    infile.getline(buffer, 90);    //! �������������Ԫ����, !before numbers there should be at least a blank
    string str1 = buffer;
    //cout << str1 << '\n';
    //sscanf_s(buffer, "%*s, %*s %d", &total_ibpoint);    //! &total_ibpoint,&ele_point);&buf, _countof(buf)
    //cout<<total_ibpoint<<endl;
    getNumberFromString(str1, a, 2);
    total_ibpoint = a[0];
    ele_point = a[1];
 
    cout << "total_ibpoint=" << total_ibpoint << "    ele_point=" << ele_point << '\n';

    lagx = new double[total_ibpoint]();
    lagy = new double[total_ibpoint]();
    lagz = new double[total_ibpoint]();

    Ep0 = new int[ele_point]();
    Ep1 = new int[ele_point]();
    Ep2 = new int[ele_point]();

    while (!infile.eof())    //! ���ж�������� ���ļ�������Ϊ��־�����ļ������һ�пհ��л�һ��δԤ�����ݣ��������ܻ�������
    {
        for (int k = 0; k < total_ibpoint; k++)
        {
            infile.getline(buffer, 100);    //! 60̫С������
            //string str1 = buffer;
            //cout << str1 << '\n';
            sscanf_s(buffer, "%lf  %lf  %lf", &lagx[k], &lagy[k], &lagz[k]);
            //cout << lagx[k] << "     " << lagy[k] << "  " << lagz[k] << '\n';
        }

        for (int k = 0; k < ele_point; k++)
        {
            infile.getline(buffer, 40);
            //string str1 = buffer;
            //cout << str1 << '\n';
            sscanf_s(buffer, "%d       %d       %d", &Ep0[k], &Ep1[k], &Ep2[k]);
            //cout << Ep0[k] << "     " << Ep1[k] << "  " << Ep2[k] << '\n';
        }
    }

    //! Ϊ�˼��һ�飬�������Finite element�������ļ�
    std::ofstream output_stream;
    output_stream.open("out.plt", std::ofstream::out);
    //! ����ʱ��option��    | std::ofstream::app���� app׷���ļ�������ΪĬ�ϸ���
    if (true)
    {
        output_stream << "Title = " << '"' << "finite-element data" << '"' << '\n';
        output_stream << "variables = x, y, z" << '\n';
        output_stream << "zone n=" << total_ibpoint << ", E=" << ele_point << ",f=FEpoint,et = triangle" << '\n';
    }

    for (int k = 0; k < total_ibpoint; k++)
    {
        output_stream << lagx[k] << "  " << lagy[k] << "  " << lagz[k] << '\n';
    }
    for (int k = 0; k < ele_point; k++)
    {
        output_stream << Ep0[k] << "  " << Ep1[k] << "  " << Ep2[k] << "  " << '\n';
    }
    output_stream.close();
 
    x_max = *max_element(lagx, lagx + total_ibpoint);
    x_min = *min_element(lagx, lagx + total_ibpoint);
    y_max = *max_element(lagy, lagy + total_ibpoint);
    y_min = *min_element(lagy, lagy + total_ibpoint);
    z_max = *max_element(lagz, lagz + total_ibpoint);
    z_min = *min_element(lagz, lagz + total_ibpoint);

    cout << "x_max=" << x_max << "  " << "xmin=" << x_min << '\n';
    cout << "y_max=" << y_max << "  " << "ymin=" << y_min << '\n';
    cout << "z_max=" << z_max << "  " << "zmin=" << z_min << '\n';
    x_length = x_max - x_min;
    y_length = y_max - y_min;
    z_length = z_max - z_min;
    std::cout <<"xlength=" << x_length << "   " << "ylength=" << y_length << "   " << "zlength=" << z_length << '\n';
    //std::cout << ac << "   " << "t=" << t << "   " << "u=" << u << "   " << "v=" << v << '\n';

    bool judge;
    vector3<double> v0, v1, v2;
    double t, u, v;
    int *t_rcd;
    t_rcd = new int[2]();
    int zmin0 = (int)floor(z_min);
    int zmax0 = (int)floor(z_max) + 2;    //! 3D
    if (direction_size == 9) { zmin0 = 0; zmax0 = 1; }     //! 2D

    for (int x = (int)floor(x_min); x < (int)floor(x_max) + 2; x++)
    {
        for (int y = (int)floor(y_min); y < (int)floor(y_max) + 2; y++)
        {
            //for (int z = (int)floor(z_min); z < (int)floor(z_max) + 2; z++) {
            for (int z =zmin0; z < zmax0; z++)
            {
                double origx = x; double origy = y; double origz = z;
                vector3<double> orig = { origx,origy,origz }, dir = { 0.0, 0.0, 1.0 };
                if (direction_size == 9) dir = { 0.0, 1.0, 0.0 };
                t_rcd[0] = 0;
                t_rcd[1] = 0;
                for (int k = 0; k < ele_point; k++)
                {
                    v0 = { lagx[Ep0[k]-1],lagy[Ep0[k]-1],lagz[Ep0[k]-1] };
                    v1 = { lagx[Ep1[k]-1],lagy[Ep1[k]-1],lagz[Ep1[k]-1] };
                    v2 = { lagx[Ep2[k]-1],lagy[Ep2[k]-1],lagz[Ep2[k]-1] };
                    if (direction_size == 9)    //! 2D
                    {
                        v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],1.0 };
                        v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],0.0 };
                        v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1], -1.0 };
                    }
                    //! ע�⣺ÿ�����ǵ�Ԫ�Ķ��㣨��㣩���Ep0[k]��Ep1[k]��Ep2[k]���Ǵ�1��ʼ�ģ��������ǵĽ����������lagx[]�Ǵ�0��ʼ�ġ�
                    //! ����ر�ע�⣺��ȷ�õ�Ep0[k]��x����ֵӦ�ò���lagx[Ep0[k]-1]
                    judge = IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &v);
                    if (judge == true )
                    {  //! && u > 0.00001 && v > 0.00001
                        if (t > 0)
                        {
                            t_rcd[0] ++;    //! tΪ��ֵ�Ĵ���
                        }
                        else
                        {
                            t_rcd[1] ++;    //! tΪ��ֵ�Ĵ��� }
                        }
                    }
                    if (t_rcd[0] > 0 && t_rcd[1] > 0) 
                    {
                        //std::cout << "p" << t_rcd[0] << "  n" << t_rcd[1] << '\n';
                        obst[index(x, y, z, NX, NY)] = 1;
                    }
                }
            }
        }
    }

    //!    ���  obst=1  tecplot�ۿ�
    std::ofstream output1_stream;
    output1_stream.open("del.plt", std::ofstream::out);
    //!    ����ʱ��option��    | std::ofstream::app���� app׷���ļ�������ΪĬ�ϸ���
    if (true) 
    {
        output1_stream << "variables = x, y, z, obst" << '\n';
        output1_stream << "zone i=" << NX << ",j=" << NY << ", k=" << NZ << ", f=point" << '\n';
    }

    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int x = 0; x < NX; x++)
            {
                output1_stream << x << " " << y << " " << z << " " <<
                    obst[index(x, y, z, NX, NY)] << '\n';
            }
        }
    }
    output1_stream.close();

    int m = 0;
    //! �������Ĺ����obstֵ�ĳ�2
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                if (obst[index(x, y, z,NX,NY)] == 1)
                {
                    int n = 0;
                    for (int k = 0; k < direction_size; k++)
                    {
                        int ind = index(x + ei[k].x, y + ei[k].y, z + ei[k].z,NX,NY);
                        if (obst[ind] == 0)
                        {
                            n++;
                        }
                    }
                    if (n > 0)
                    {
                        obst[index(x, y, z, NX, NY)] = 2;
                        m++;
                    }
                }
            }
        }
    }   //! Ŀ���ǵõ�obst2�Ľ����  ��ȷ��������Ĵ�С

    num_node2 = m;

    m = 0;
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                if (obst[index(x, y, z,NX,NY)] == 2)
                {
                    node2[m] = index(x, y, z, NX, NY);
                    for (int k1 = 0; k1 < direction_size; k1++)
                    {
                        int ip = ei[k1].x, jp = ei[k1].y, kp = ei[k1].z;
                        int ind = index(x + ip, y + jp, z + kp, NX, NY);
                        bool judge;
                        double t, u, v;
                        vector3<double> v0, v1, v2;
                        vector3<double> orig = { (double)x,(double)y,(double)z },
                            dir = { (double)ip, (double)jp, (double)kp };
                        double b1 = 0., a1 =  dir.norm_square();    //! �ٶ�ģ�����ٶȳ��� 1, sqrt2, sqrt3...
                        if (obst[ind] == 0)
                        {
                            for (int k = 0; k < ele_point; k++)
                            {
                                v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],lagz[Ep0[k] - 1] };
                                v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],lagz[Ep1[k] - 1] };
                                v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1],lagz[Ep2[k] - 1] };
                                if (direction_size == 9)
                                {
                                    v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],1.0 };
                                    v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],0.0 };
                                    v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1],-1.0 };
                                }
                                judge = IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &v);
                                if (judge == true && t>0 )
                                {
                                    b1 = t;
                                }    //! && t<b1
                            }
                            if (b1 / a1 > 1) delta[m * direction_size + k1] = 0.;
                            else delta[m * direction_size + k1] = 1.0- b1/sqrt(a1);
                            //cout <<"q="<< 1.0- b1/a1<<"    " <<"a1=" <<a1 <<'\n';
                            //if (b1 / a1 > 1) cout <<"b1="<<b1<<" a1="<<a1<<" error!!!!!!"<< '\n';
                        }
                    }
                    m++;
                }
            }
        }
    }
    //int aaa;
}



void getNumberFromString(string s, int *a, int n) 
{
    stringstream str_strm;
    str_strm << s;    //! convert the string s into stringstream
    string temp_str;
    int temp_int;
    int k = 0;
    while (!str_strm.eof())
    {
        str_strm >> temp_str;    //! take words into temp_str one by one
        if (stringstream(temp_str) >> temp_int)    //! try to convert string to int
        {
            a[k] = temp_int;    cout << temp_int << endl;
            k++;
        }
        temp_str = "";    //! clear temp string
    }
}


bool IntersectTriangle( vector3<double>& orig, vector3<double>& dir,
    vector3<double>& v0, vector3<double>& v1, vector3<double>& v2,
    double* t, double* u, double* v)
{
    //! E1
    vector3<double> E1 = v1 - v0;

    //! E2
    vector3<double> E2 = v2 - v0;

    //! P
    vector3<double> P = dir.cross(E2);

    //! determinant
    double det = E1.dot(P);

    //! keep det > 0, modify T accordingly
    vector3<double> T;
    if (det > 0)
    {
        T = orig - v0;
    }
    else
    {
        T = v0 - orig;
        det = -det;
    }

    //! If determinant is near zero, ray lies in plane of triangle
    if (det < 0.0001f)
        return false;

    //! Calculate u and make sure u <= 1
    *u = T.dot(P);
    if (*u < 0.0f || *u > det)
        return false;

    //! Q
    vector3<double> Q = T.cross(E1);

    //! Calculate v and make sure u + v <= 1
    *v = dir.dot(Q);
    if (*v < 0.0f || *u + *v > det)
        return false;

    //! Calculate t, scale parameters, ray intersects triangle
    *t = E2.dot(Q);

    double fInvDet = 1.0f / det;
    *t *= fInvDet;
    *u *= fInvDet;
    *v *= fInvDet;

    return true;
}