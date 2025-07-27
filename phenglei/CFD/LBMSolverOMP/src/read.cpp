#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "LBMSolverOMP.hpp"
using namespace std;
//******************************************
// R E A D   P A R A M S (simulation parameters)
//
//  - Read the problem parameters from a file.
//

void read_params(Params_Pakage* p)
{
    char buffer[256];
    std::string A;    //! �ַ���
    vector3<double> aa;
    double rho_specified;
    ifstream infile(p->in);     //! �����ļ�
    ofstream outfile(p->out);   //! ����ļ�

    if (!infile)
    {
        cout << "Unable to open input file";
        exit(1); // terminate with error
    }
    if (!outfile)
    {
        cout << "Unable to open output file";
        exit(1); // terminate with error
    }
    //sscanf_s(tokenstring, "%s", s, _countof(s));
    //sscanf_s(tokenstring, "%c", &c, sizeof(char));

    while (!infile.eof())    //! ���ж�������� ���ļ�������Ϊ��־�����ļ������һ�пհ��л�һ��δԤ�����ݣ��������ܻ�������
    {
        infile.getline(buffer, 100);    //! �����־�Ƿ����㣬��Ϊ1����Ϊ0������������������������ļ���
        sscanf_s(buffer, "%d", &p->MRT );
        infile.getline(buffer, 120);    //! �����־�Ƿ����㣬��Ϊ1����Ϊ0������������������������ļ���
        sscanf_s(buffer, "%d", &p->CONTI);
        if(p->CONTI==1) sscanf_s(buffer, "%s", &p->conti_file, (unsigned int)sizeof(p->conti_file));
        infile.getline(buffer, 100);    //! �����־�Ƿ����ģ�⣬��Ϊ1����Ϊ0
        sscanf_s(buffer, "%d", &p->LES);
        infile.getline(buffer, 120);    //! �����־�Ƿ���Ҫ���븴���������壬��Ϊ1����Ϊ0����Ҫ���������ͼ��tecplot�����ļ���
        sscanf_s(buffer, "%d", &p->GEO);
        if (p->GEO == 1) sscanf_s(buffer, "%s", &p->geo_file, (unsigned int)sizeof(p->geo_file));
        infile.getline(buffer, 100);    //! �����־�Ƿ���������Ϊ1����Ϊ0
        sscanf_s(buffer, "%d", &p->MB);
        infile.getline(buffer, 100);    //! �����־�Ƿ���Ҫ�����оݣ���Ϊ1����Ϊ0������Ҫ���������������ֵ
        sscanf_s(buffer, "%d", &p->CONV);
        if (p->CONV == 1) sscanf_s(buffer, "%lf", &p->conv_v);

        infile.getline(buffer, 100);    //! ��������� X ����������
        sscanf_s(buffer, "%d", &p->LX);
        infile.getline(buffer, 100);    //! ������ Y ����������
        sscanf_s(buffer, "%d", &p->LY);
        infile.getline(buffer, 100);    //! ������ Z ����������
        sscanf_s(buffer, "%d", &p->LZ);
        infile.getline(buffer, 100);    //! ������ X ����������
        sscanf_s(buffer, "%d", &p->LX2);
        infile.getline(buffer, 100);    //! ������ Y ����������
        sscanf_s(buffer, "%d", &p->LY2);
        infile.getline(buffer, 100);    //! ������ Z ����������
        sscanf_s(buffer, "%d", &p->LZ2);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowX);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowY);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowZ);
        infile.getline(buffer, 100);    //! ÿ�����ٲ����һ���ļ�

        sscanf_s(buffer, "%d", &p->FrameRate);
        infile.getline(buffer, 100);    //! �ܹ����ٲ���ֹ����
        sscanf_s(buffer, "%d", &p->ttstep);
        infile.getline(buffer, 100);    //! ��ʼ�ܶ�
        sscanf_s(buffer, "%lf", &p->density0);
        infile.getline(buffer, 100);    //! �ɳ����ӣ����˶�ճ�Թ�ϵ:\nu=c_s^2*(\tau-0.5)\delta t
        sscanf_s(buffer, "%lf", &p->tau);
        infile.getline(buffer, 120);    //! �����ʸ��
        sscanf_s(buffer, "(%lf, %lf, %lf)", &p->gf.x, &p->gf.y, &p->gf.z);
        infile.getline(buffer, 100);    //! LBM�ٶ�ģ��
        sscanf_s(buffer, "%s", &p->velocity_set, (unsigned int)sizeof(p->velocity_set));
        //! �ַ�����ȡʱ��Ĭ���Կո�ָ�
        infile.getline(buffer, 100);    //! �߽���������
        sscanf_s(buffer, "%s", &p->bc, (unsigned int)sizeof(p->bc));
        infile.getline(buffer, 256);    //! �ļ��е�˵���У���ȡ��Ϣ��������
        infile.getline(buffer, 256);    //! �ļ��е�˵���У���ȡ��Ϣ��������

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.xminFace, (unsigned int)sizeof(p->BC.xminFace));
        A = p->BC.xminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[0] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[0] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.xmaxFace, (unsigned int)sizeof(p->BC.xmaxFace));
        A = p->BC.xmaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[1] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[1] = rho_specified;
        }
        
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.yminFace, (unsigned int)sizeof(p->BC.yminFace));
        A = p->BC.yminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[2] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[2] = rho_specified;
        }

        infile.getline(buffer, 256);
        sscanf_s(buffer, "%s", &p->BC.ymaxFace, (unsigned int)sizeof(p->BC.ymaxFace));
        A = p->BC.ymaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[3] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[3] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.zminFace, (unsigned int)sizeof(p->BC.zminFace));
        A = p->BC.zminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[4] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[4] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.zmaxFace, (unsigned int)sizeof(p->BC.zmaxFace));
        A = p->BC.zmaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[5] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[5] = rho_specified;
        }

        infile.getline(buffer, 100);    //! ��һ��
        infile.getline(buffer, 100);
        sscanf_s(buffer, "(%lf, %lf, %lf)", &p->vel0.x, &p->vel0.y, &p->vel0.z);

        infile.getline(buffer, 100);    //!IBFSIģ���Ƿ�ʹ�ö���ṹ
        sscanf_s(buffer, "%d", &p->IB_Multi_Geo);
        infile.getline(buffer, 100);    //! �Ƿ�ʹ��IBFSIģ��
        sscanf_s(buffer, "%d", &p->IBFSI);

        cout <<"Computational domain size:   "<< p->LX  <<"  "<< p->LY<<"   "<< p->LZ<< endl;
        cout <<"The velocity model used:   "<< p->velocity_set <<'\n'<< endl;

        cout << p->bc << '\n' << endl;

    }

    infile.close();

    //! ��ͬʱ������ϲ�������һ�ļ�������֤��ȡ����ȷ��
    outfile << "MRT:" << "   " << p->MRT << endl;
    outfile << "CONTI:" << "   " << p->CONTI << endl;
    outfile << "LES:" << "   " << p->LES << endl; 
    outfile << "GEO:" << "   " << p->GEO << endl;
    outfile << "MB:" << "   " << p->MB << endl;
    outfile << "LX:" << "   " << p->LX << endl;
    outfile << "LY:" << "   " << p->LY << endl;
    outfile << "LZ:" << "   " << p->LZ << endl;
    outfile << "LX2:" << "   " << p->LX2 << endl;
    outfile << "LY2:" << "   " << p->LY2<< endl;
    outfile << "LZ2:" << "   " << p->LZ2 << endl;
    outfile << "LowX:" << "   " << p->LowX << endl;
    outfile << "LowY:" << "   " << p->LowY << endl;
    outfile << "LowZ:" << "   " << p->LowZ << endl;
    outfile << "FrameRate:" << "   " << p->FrameRate<< endl;
    outfile << "Total steps:" << " " << p->ttstep << endl;
    outfile << "LBM.initial density:" << " " << p->density0 << endl;
    if (p->IBFSI) 
    {
        outfile << "LBM.tau:" << " Calculated by Re in IBFSI" << endl;
    }
    else 
    {
        outfile << "LBM.tau:" << " " << p->tau << endl;
    }
    outfile << "LBM.tau:" << " " << p->tau << endl;
    outfile << "LBM.gx:" << " " << p->gf.x << endl;
    outfile << "LBM.gy:" << " " << p->gf.y << endl;
    outfile << "LBM.gz:" << " " << p->gf.z << endl;
    outfile << "velocity model:" << " " << p->velocity_set << endl;
    outfile << "Boundary condition type:" << " " << p->bc << endl;
    outfile << "xmin Face BC:  " << " " << p->BC.xminFace << endl;
    outfile << "xmax Face BC:  " << " " << p->BC.xmaxFace << endl;
    outfile << "ymin Face BC:  " << " " << p->BC.yminFace << endl;
    outfile << "ymax Face BC:  " << " " << p->BC.ymaxFace << endl;
    outfile << "zmin Face BC:  " << " " << p->BC.zminFace << endl;
    outfile << "zmax Face BC:  " << " " << p->BC.zmaxFace << endl;
    outfile << "initial vel  " << " " << p->vel0.x << " " << p->vel0.y << " " << p->vel0.z << endl;
    outfile << "multiple geo number:" << " " << p->IB_Multi_Geo << endl;
    outfile << "UseIBPeskin or not:" << " " << p->IBFSI << endl;

    outfile.close();

//    return 0;
}





