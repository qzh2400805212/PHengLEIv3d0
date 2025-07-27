#include "LBMSolverOMP.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <iostream>


void LBM::wall_bounce_back()
{   //! nonslip���棬bounce-backʵ���޻��Ʊ߽�����.
    std::string ymin = BC.yminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ���׼��
    if (ymin == "nonslip")
    {
        int y = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                   if (ei[i].y == 1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }
    std::string ymax = BC.ymaxFace;
    if (ymax == "nonslip")
    {
        int y = NY-1;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                 if (ei[i].y == -1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }
    //! x direction
    std::string xmin = BC.xminFace;
    if (xmin == "nonslip")
    {
        int x = 0;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    if (ei[i].x == 1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }
    std::string xmax = BC.xmaxFace;
    if (xmax == "nonslip")
    {
        int x = NX - 1;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    if (ei[i].x == -1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }
    //! z direction
    std::string zmin = BC.zminFace;
    if (zmin == "nonslip")
    {
        int z = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    if (ei[i].z == 1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }
    std::string zmax = BC.zmaxFace;
    if (zmax == "nonslip")
    {
        int z = NZ - 1;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    if (ei[i].z == -1) pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                }
            }
        }
    }

}

void LBM::velocity_boundary_condition()
{
    const double tauinv = 1.0 / tau;
    const double omtauinv = 1.0 - tauinv;

    std::string Astring = BC.xminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "velocity")
    {
        int x = 0;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[0].x;
                vel[index(x, y, z)].y = BC.vel_specified[0].y;
                vel[index(x, y, z)].z = BC.vel_specified[0].z;
                rho[index(x, y, z)] = rho[index(x+1, y, z)];    //! rho �����������������(x+1, y, z)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x+1, y, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x+1, y, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x+1, y, z, i)] - feq1) �����������������(x+1, y, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.xmaxFace;
    if (Astring == "velocity")
    {
        int x = NX-1;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[1].x;
                vel[index(x, y, z)].y = BC.vel_specified[1].y;
                vel[index(x, y, z)].z = BC.vel_specified[1].z;
                rho[index(x, y, z)] = rho[index(x - 1, y, z)];    //! rho �����������������(x-1, y, z)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x - 1, y, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x - 1, y, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x-1, y, z, i)] - feq1) �����������������(x-1, y, z)��ֵ����
                }
            }
        }
    }



    Astring = BC.yminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    /*if (Astring == "velocity"&& BC.vel_specified[2].y == 0.0) {//�������棺�����ٶ� 
        int y = 0;
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int i = 0; i < direction_size; i++) {
                    if (ei[i].y == 1) {
                        pdf[index(x, y, z, i)] = pdf2[index(x, y, z, reverse_indexes[i])] 
                            + 2 * weights[i] / (c_s * c_s) *( ei[i].x * BC.vel_specified[2].x 
                                + ei[i].z * BC.vel_specified[2].z);
                    }
                }
            }
        }
    }*/


    //! ���¸���������ϵ��ٶȣ����������ⷽ�򣨰�����ֱ�ڸ���Ľ����ٶȣ�Ҳ�����ǻ������棨ƽ�и�����ٶȣ���
    if (Astring == "velocity")
    {
        int y = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[2].x;
                vel[index(x, y, z)].y = BC.vel_specified[2].y;
                vel[index(x, y, z)].z = BC.vel_specified[2].z;
                rho[index(x, y, z)] = rho[index(x, y + 1, z)];    //! rho �����������������(x, y+1, z)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y + 1, z, i);
                    pdf[index(x, y, z, i)] =feq +omtauinv * (pdf2[index(x, y+1, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y+1, z, i)] - feq1) �����������������(x, y+1, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.ymaxFace;
/*    if (Astring == "velocity"&& BC.vel_specified[3].y == 0.0) {
        int y = NY - 1;
        for (int x = 0; x < NX; x++) {
            for (int z = 0; z < NZ; z++) {
                for (int i = 0; i < direction_size; i++) {//�������棺�����ٶ� 
                    if (ei[i].y == -1) {
                        
                        pdf[index(x, y, z, i)] = pdf2[index(x, y, z, reverse_indexes[i])] 
                            + 2 * weights[i] / (c_s * c_s) *(ei[i].x * BC.vel_specified[3].x
                                + ei[i].z * BC.vel_specified[3].z);
                    }
                }
            }
        }
    }*/
    if (Astring == "velocity")
    {
        int y = NY-1;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[3].x;
                vel[index(x, y, z)].y = BC.vel_specified[3].y;
                vel[index(x, y, z)].z = BC.vel_specified[3].z;
                rho[index(x, y, z)] = rho[index(x, y - 1, z)];    //! rho �����������������(x, y-1, z)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y - 1, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y - 1, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y+1, z, i)] - feq1) �����������������(x, y+1, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.zminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "velocity")
    {
        int z = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[4].x;
                vel[index(x, y, z)].y = BC.vel_specified[4].y;
                vel[index(x, y, z)].z = BC.vel_specified[4].z;
                rho[index(x, y, z)] = rho[index(x, y, z+1)];    //! rho �����������������(x, y, z+1)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y, z+1, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z+1, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y, z+1, i)] - feq1) �����������������(x, y, z+1)��ֵ����
                }
            }
        }
    }

    Astring = BC.zmaxFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "velocity")
    {
        int z = NZ-1;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                vel[index(x, y, z)].x = BC.vel_specified[5].x;
                vel[index(x, y, z)].y = BC.vel_specified[5].y;
                vel[index(x, y, z)].z = BC.vel_specified[5].z;
                rho[index(x, y, z)] = rho[index(x, y, z - 1)];    //! rho �����������������(x, y, z-1)��ֵ����
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y, z-1, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z - 1, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y, z-1, i)] - feq1) �����������������(x, y, z-1)��ֵ����
                }
            }
        }
    }

}



void LBM::pressure_boundary_condition()
{
    const double tauinv = 1.0 / tau;
    const double omtauinv = 1.0 - tauinv;

    std::string Astring = BC.xminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    //! ���¸������ϵ�ѹǿ���ܶȣ�
    if (Astring == "pressure")
    {
        int x = 0;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = vel[index(x + 1, y, z)].x;    //! �ٶ������������������(x+1, y, z)��ֵ����
                vel[index(x, y, z)].y = vel[index(x + 1, y, z)].y;
                vel[index(x, y, z)].z = vel[index(x + 1, y, z)].z;
                rho[index(x, y, z)] = BC.rho_specified[0];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x + 1, y, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x + 1, y, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x+1, y, z, i)] - feq1) �����������������(x+1, y, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.xmaxFace;
    if (Astring == "pressure")
    {
        int x = NX - 1;
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = vel[index(x - 1, y, z)].x;
                vel[index(x, y, z)].y = vel[index(x - 1, y, z)].y;
                vel[index(x, y, z)].z = vel[index(x - 1, y, z)].z;
                rho[index(x, y, z)] = BC.rho_specified[1];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x - 1, y, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x - 1, y, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x-1, y, z, i)] - feq1) �����������������(x-1, y, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.yminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "pressure")
    {
        int y = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = vel[index(x, y+1, z)].x;
                vel[index(x, y, z)].y = vel[index(x, y+1, z)].y;
                vel[index(x, y, z)].z = vel[index(x, y+1, z)].z;
                rho[index(x, y, z)] = BC.rho_specified[2];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y + 1, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y + 1, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y+1, z, i)] - feq1) �����������������(x, y+1, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.ymaxFace;
    if (Astring == "pressure")
    {
        int y = NY - 1;
        for (int x = 0; x < NX; x++)
        {
            for (int z = 0; z < NZ; z++)
            {
                vel[index(x, y, z)].x = vel[index(x, y - 1, z)].x;
                vel[index(x, y, z)].y = vel[index(x, y - 1, z)].y;
                vel[index(x, y, z)].z = vel[index(x, y - 1, z)].z;
                rho[index(x, y, z)] = BC.rho_specified[3];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y - 1, z, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y - 1, z, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y+1, z, i)] - feq1) �����������������(x, y+1, z)��ֵ����
                }
            }
        }
    }

    Astring = BC.zminFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "pressure")
    {
        int z = 0;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                vel[index(x, y, z)].x = vel[index(x, y, z+1)].x;
                vel[index(x, y, z)].y = vel[index(x, y, z+1)].y;
                vel[index(x, y, z)].z = vel[index(x, y, z+1)].z;
                rho[index(x, y, z)] = BC.rho_specified[4];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y, z + 1, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z + 1, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y, z+1, i)] - feq1) �����������������(x, y, z+1)��ֵ����
                }
            }
        }
    }

    Astring = BC.zmaxFace;    //! �ַ����鸳ֵ��string��Ϊ�ַ����Ƚ�׼��
    if (Astring == "pressure")
    {
        int z = NZ - 1;
        for (int x = 0; x < NX; x++)
        {
            for (int y = 0; y < NY; y++)
            {
                vel[index(x, y, z)].x = vel[index(x, y, z - 1)].x;    //! �����������������(x, y, z - 1)��ֵ����
                vel[index(x, y, z)].y = vel[index(x, y, z - 1)].y;
                vel[index(x, y, z)].z = vel[index(x, y, z - 1)].z;
                rho[index(x, y, z)] = BC.rho_specified[5];
                for (int i = 0; i < direction_size; i++)
                {
                    double feq = calculate_feq(x, y, z, i);
                    double feq1 = calculate_feq(x, y, z - 1, i);
                    pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z - 1, i)] - feq1);
                    //! ��ƽ�ⲿ��(pdf2[index(x, y, z-1, i)] - feq1) �����������������(x, y, z-1)��ֵ����
                }
            }
        }
    }

}
