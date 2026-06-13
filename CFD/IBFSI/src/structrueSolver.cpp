#include "IB.hpp"

void trid(int N, vector<double> a, vector<double> b, vector<double> c, vector<double>& x, vector<double> f);
void pdma2(int n, vector<double> fAww, vector<double> fAw, vector<double> fAp, vector<double> fAe,
    vector<double> fAee, vector<double>& x, vector<double> frhs);


void IBFSI::prescribed_motion(int time, int type)
{
    if (type == 1)    //! flow past a fixed cylinder
    {
        Xnode = X0node;
        Ynode = Y0node;
        Unode = creatmatrix_1D(Nnode, 0.0);
        Vnode = creatmatrix_1D(Nnode, 0.0);
    }
    else if (type == 2)    //! flapping rigid plate
    {
        for (int i = 0; i < Nnode; ++i)
        {
            int rotateOrnot = IsRotateOrnot;
            double theta0 = iniAngle * pi / 180.0, theta0_ = iniPhase * pi / 180.0;    //! 初始角度与角度的初始相位
            double A_theta = A_rotate;    //! 刚性板转动的赋值：尾缘的幅度
            double T_theta = T_rotate;    //! 拍动频率：与特征周期的关系

            double theta_max = asin(A_theta / Lref);    //! 使用尾缘幅度求出相应弧度
            double TT = Lref / (ds0 * Uref);    //! 特征周期
            double f = 1 / (T_theta * TT);
            double theta = theta0 + theta_max * cos(2 * pi * f * time + theta0_);
            double sintheta = sin(2 * pi * f * time + theta0_);

            if (rotateOrnot == 0)
            {
                f = 0;
                theta = theta0;
            }

            Xnode[i] = Xnode[0] + i * ds0 * cos(theta);
            Ynode[i] = Ynode[0] + i * ds0 * sin(theta);
            Unode[i] = i * ds0 * 2 * pi * A_theta * f * sin(theta) * sintheta;
            Vnode[i] = -1 * (i * ds0 * 2 * pi * A_theta * f * cos(theta) * sintheta);
        }
    }
    else
    {
        cout << "error prescribed_motion: please input the type 1 or 2!!!" << endl;
    }
}


void IBFSI::filament_motion(int iternum)
{
    vector<double> Fr = { gravityX,gravityY,gravityZ };    //! fr = gL / U ^ 2 注意正负
    //vector<double> Fr = { 0.0,-50,0.0 };
 
    //! 强制运动控制参数
    double A_heaving = Ampl_heaving;    //! 强制位移幅值
    double T_heaving = Period_heaving;    //! 强制运动周期

    double A_pitching = Ampl_pitching;    //! 强制转角幅值
    double T_pitching = Period_pitching;    //! 强制运动周期
    double theta_ini = asin(A_pitching / Lref);

    //! heaving
    vector<double> velfix(2, 0.0);
    double f_heaving = (T_heaving != 0) ? (1.0 / (T_heaving * Tref)) : 0;
    double f_pitching = (T_pitching != 0) ? 1.0 / (T_pitching * Tref) : 0;
    velfix[0] = 0.0;
    velfix[1] = 2.0 * pi * A_heaving * f_heaving * sin(2.0 * pi * f_heaving * solvertime);

    //! pitching
    vector<double> angfix(2, 0.0);
    double theta = theta_ini * cos(2 * pi * f_pitching * solvertime);
    angfix[0] = cos(theta);
    angfix[1] = sin(theta);

    vector<vector<double>> xfo(Nnode, vector<double>(2, 0.0));
    vector<vector<double>> xfn(Nnode, vector<double>(2, 0.0));
    vector<vector<double>> ufn(Nnode, vector<double>(2, 0.0));
    vector<vector<double>> exf(Nnode, vector<double>(2, 0.0));

    if (solvertime == 0)
    {
        for (int i = 0; i < Nnode; i++)
        {
            Xnode[i] = Xnode[0] + i * ds0 * cos(theta);
            Ynode[i] = Ynode[0] + i * ds0 * sin(theta);
        }
        X0node = Xnode;
        Y0node = Ynode;
    }

    for (int ntemp = 0; ntemp < Nnode; ++ntemp)
    {
        xfo[ntemp][0] = X0node[ntemp] / Lref;
        xfo[ntemp][1] = Y0node[ntemp] / Lref;
        xfn[ntemp][0] = Xnode[ntemp] / Lref;
        xfn[ntemp][1] = Ynode[ntemp] / Lref;
        ufn[ntemp][0] = Unode[ntemp] / Uref;
        ufn[ntemp][1] = Vnode[ntemp] / Uref;
        exf[ntemp][0] = (extfulx[ntemp] + Fr[0]) / Fref;
        exf[ntemp][1] = (extfuly[ntemp] + Fr[1]) / Fref;
        //exf[ntemp][0] = Fr[0];
        //exf[ntemp][1] = Fr[1];
    }

    for (int i = 0; i < int(iternum * Tref); ++i)
    {
        filament_solver(dt/(Tref * double(iternum)), velfix, angfix, xfo, xfn, ufn, exf);
    }

    solvertime += dt;

    for (int ntemp = 0; ntemp < Nnode; ++ntemp)
    {
        X0node[ntemp] = xfo[ntemp][0] * Lref;
        Y0node[ntemp] = xfo[ntemp][1] * Lref;
        Xnode[ntemp] = xfn[ntemp][0] * Lref;
        Ynode[ntemp] = xfn[ntemp][1] * Lref;
        Unode[ntemp] = ufn[ntemp][0] * Uref;
        Vnode[ntemp] = ufn[ntemp][1] * Uref;
    }
}

void IBFSI::VIV_motion()
{
    double CylinderForce_y = 0.0;
    for (int i = 0; i < Nnode; ++i)
    {
        CylinderForce_y += Fynode[i];
    }
    VIVForce_Fy = CylinderForce_y;

    VIV_Qn[0] = VIVForce_Fy;
    VIV_Qn[1] = 0.0;    //! 无力矩

    solve_VIV(dt, VIV_X, VIV_Fn, VIV_Fn_1, VIV_Fn_2, VIV_Fn_3,
        VIV_Qn, VIV_Qn_1, VIV_Qn_2, VIV_Qn_3, VIV_A, VIV_B);

    for (int n = 0; n < Nnode; ++n)
    {
        Xnode[n] = VIVMassCenter_X 
                 + ((X00node[n] - VIVMassCenter_X) * cos(VIV_X[1]) 
                  - (Y00node[n] - VIVMassCenter_Y) * sin(VIV_X[1]));
        Ynode[n] = VIVMassCenter_Y + VIV_X[0] 
                 + ((X00node[n] - VIVMassCenter_X) * sin(VIV_X[1]) 
                  + (Y00node[n] - VIVMassCenter_Y) * cos(VIV_X[1]));

        Unode[n] = 0.0 - VIV_X[3] * (Ynode[n] - VIVMassCenter_Y);
        Vnode[n] = VIV_X[2] + VIV_X[3] * (Xnode[n] - VIVMassCenter_X);
        Wnode[n] = 0.0;
    }
}

void IBFSI::filament_solver(double dtfilm, vector<double> velfix, vector<double> angfix, vector<vector<double>> &xfo, vector<vector<double>> &xfn,
    vector<vector<double>> &ufn, vector<vector<double>> exf)
{
    int nfilm = Nelem;
    double dsfilm = ds0;

    int isFree = 1;

    int ifixend = subMotionType1;
    int iextend = subMotionType2;
    double kbend = kb_film;
    double kstrech = ks_film;
    //double kbend = kb_film / KBref;
    //double kstrech = ks_film / KSref;
    double mass_beta = 1.0;

    //vector<double> angfix = { angfix_film1,angfix_film2 };
    vector<double> accfix = { accfix_film1,accfix_film2 };

    //! Coefficients and rhs term for Equation (26) in JCP paper by Huang WeiXi
    vector<vector<double>> fAww(nfilm + 1, vector<double>(2, 0.0));   //! Coefficient matrix Aww
    vector<vector<double>> fAw(nfilm + 1, vector<double>(2, 0.0));    //! Coefficient matrix Aw
    vector<vector<double>> fAp(nfilm + 1, vector<double>(2, 0.0));    //! Coefficient matrix Ap
    vector<vector<double>> fAe(nfilm + 1, vector<double>(2, 0.0));    //! Coefficient matrix Ae
    vector<vector<double>> fAee(nfilm + 1, vector<double>(2, 0.0));   //! Coefficient matrix Aee
    vector<vector<double>> bfilm(nfilm + 1, vector<double>(2, 0.0));  //! Right-hand side vector b

    vector<double> TAw(nfilm + 1, 0.0);    //! Coefficients TAw
    vector<double> TAp(nfilm + 1, 0.0);    //! Coefficients TAp
    vector<double> TAe(nfilm + 1, 0.0);    //! Coefficients TAe
    vector<double> rhsT(nfilm, 0.0);       //! Right-hand side vector rhsT

    //! Lagrangian coordinates of the nodes on the filament
    vector<vector<double>> xft(nfilm + 1, vector<double>(2, 0.0));    //! Lagrangian coordinates xft
    vector<vector<double>> xfa(nfilm + 1, vector<double>(2, 0.0));    //! Lagrangian coordinates xfa (x*)

    //! First and second derivative of node position
    vector<vector<double>> dx1s(nfilm, vector<double>(2, 0.0));    //! First derivative dx1s
    vector<vector<double>> dx2s(nfilm + 1, vector<double>(2, 0.0));    //! Second derivative dx2s

    //! Tension, bending force, and first derivative of exerting force
    vector<double> tension(nfilm, 0.0);    //! Tension
    vector<vector<double>> bendf(nfilm + 1, vector<double>(2, 0.0));    //! Bending force bendf
    vector<vector<double>> df1s(nfilm, vector<double>(2, 0.0));    //! First derivative of exerting force df1s

    double eps, teps;     //! Variables eps and teps

    //! Constants for calculations
    double pdsfilm = 1.0 / dsfilm;    //! Inverse of dsfilm
    double p2dsfilm = 1.0 / (dsfilm * dsfilm);    //! Inverse of dsfilm squared
    double p3dsfilm = pdsfilm * p2dsfilm;    //! Inverse of dsfilm cubed
    double p4dsfilm = p2dsfilm * p2dsfilm;    //! Fourth power of the inverse of dsfilm
    double pdtfilm = 1.0 / dtfilm;    //! Inverse of dtfilm
    double p2dtfilm = 1.0 / (dtfilm * dtfilm);    //! Inverse of dtfilm squared

    //! 确定不同端点处节点的加速度二阶导数和固定端点处的条件
    //! Calculate xft: Lagrangian coordinates of the nodes on the filament
    //! xft 预测的下一时刻的拉格朗日位移
    for (int ns = 0; ns <= nfilm; ++ns)
    {
        xft[ns][0] = 2.0 * xfn[ns][0] - xfo[ns][0];
        xft[ns][1] = 2.0 * xfn[ns][1] - xfo[ns][1];
    }

    //! Calculate position of the fixed end using velocity and acceleration boundary conditions
    xft[0][0] = xfn[0][0] + velfix[0] * dtfilm + 0.5 * accfix[0] * dtfilm * dtfilm;
    xft[0][1] = xfn[0][1] + velfix[1] * dtfilm + 0.5 * accfix[1] * dtfilm * dtfilm;

    // Calculate the second derivative of X at i
    if (ifixend == 1)    //! Simply supported: d^2X/ds^2 = (0,0)
    {
        dx2s[0][0] = 0.0;
        dx2s[0][1] = 0.0;
    }
    else if (ifixend == 2)    //! Clamped: dX/ds = (cos, sin)
    {
        dx2s[0][0] = ((xft[1][0] - xft[0][0]) * pdsfilm - angfix[0]) * 2.0 * pdsfilm;
        dx2s[0][1] = ((xft[1][1] - xft[0][1]) * pdsfilm - angfix[1]) * 2.0 * pdsfilm;
    }

    //! Calculate the second derivative of X at intermediate nodes
    for (int ns = 1; ns < nfilm; ++ns)
    {
        dx2s[ns][0] = (xft[ns + 1][0] - 2.0 * xft[ns][0] + xft[ns - 1][0]) * p2dsfilm;
        dx2s[ns][1] = (xft[ns + 1][1] - 2.0 * xft[ns][1] + xft[ns - 1][1]) * p2dsfilm;
    }

    if (isFree == 1)    //! Free end: d^2X/ds^2 = (0,0)
    {
        dx2s[nfilm][0] = 0.0;
        dx2s[nfilm][1] = 0.0;
    }
    else if (isFree == 2)    //! Simply supported: d^2X/ds^2 = (0,0)
    {
        dx2s[nfilm][0] = 0.0;
        dx2s[nfilm][1] = 0.0;
    }

    //! 计算弯曲力
    //! Calculate bending force at the first node
    bendf[0][0] = -kbend * (-dx2s[3][0] + 4.0 * dx2s[2][0] - 5.0 * dx2s[1][0] + 2.0 * dx2s[0][0]) * p2dsfilm;
    bendf[0][1] = -kbend * (-dx2s[3][1] + 4.0 * dx2s[2][1] - 5.0 * dx2s[1][1] + 2.0 * dx2s[0][1]) * p2dsfilm;

    //! Calculate bending force at intermediate nodes
    for (int ns = 1; ns < nfilm; ++ns) {
        bendf[ns][0] = -kbend * (dx2s[ns + 1][0] - 2.0 * dx2s[ns][0] + dx2s[ns - 1][0]) * p2dsfilm;
        bendf[ns][1] = -kbend * (dx2s[ns + 1][1] - 2.0 * dx2s[ns][1] + dx2s[ns - 1][1]) * p2dsfilm;
    }

    //! Calculate bending force at the last node
    //bendf[nfilm][0] = -kbend * (0.0 - (dx2s[nfilm - 1][0] - dx2s[nfilm - 2][0]) * pdsfilm) * pdsfilm; // Free end: d^3X/ds^3 = (0,0)
    //bendf[nfilm][1] = -kbend * (0.0 - (dx2s[nfilm - 1][1] - dx2s[nfilm - 2][1]) * pdsfilm) * pdsfilm; // Free end: d^3X/ds^3 = (0,0)
    if (isFree == 1)    //! Free end: d^2X/ds^2 = (0,0)
    {
        bendf[nfilm][0] = -kbend * (0.0 - (dx2s[nfilm - 1][0] - dx2s[nfilm - 2][0]) * pdsfilm) * pdsfilm; // Free end: d^3X/ds^3 = (0,0)
        bendf[nfilm][1] = -kbend * (0.0 - (dx2s[nfilm - 1][1] - dx2s[nfilm - 2][1]) * pdsfilm) * pdsfilm; // Free end: d^3X/ds^3 = (0,0)
    }
    else if (isFree == 2)    //! fixed supported
    {
        bendf[nfilm][0] = kbend * (-dx2s[nfilm-3][0] + 4.0 * dx2s[nfilm-2][0] - 5.0 * dx2s[nfilm-1][0] + 2.0 * dx2s[nfilm][0]) * p2dsfilm;
        bendf[nfilm][1] = -kbend * (-dx2s[nfilm-3][1] + 4.0 * dx2s[nfilm-2][1] - 5.0 * dx2s[nfilm-1][1] + 2.0 * dx2s[nfilm][1]) * p2dsfilm;
    }

    //! 计算中间节点拉伸力
    //! Calculate tension force at intermediate time step
    if (iextend == 1)
    {
        //! Solve constitutive equation for T
        for (int ns = 0; ns < nfilm; ++ns)
        {
            dx1s[ns][0] = (xft[ns + 1][0] - xft[ns][0]) * pdsfilm;
            dx1s[ns][1] = (xft[ns + 1][1] - xft[ns][1]) * pdsfilm;

            tension[ns] = kstrech * (sqrt(dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1]) - 1.0);
        }
        tension[nfilm] = 0.0;
    }
    else if (iextend == 2)
    {
        //! Solve Poisson equation for T
        for (int ns = 0; ns < nfilm; ++ns)
        {   //! (df1s * dx1s) at i+1/2
            df1s[ns][0] = (bendf[ns + 1][0] + exf[ns + 1][0] - bendf[ns][0] - exf[ns][0]) * pdsfilm;
            df1s[ns][1] = (bendf[ns + 1][1] + exf[ns + 1][1] - bendf[ns][1] - exf[ns][1]) * pdsfilm;

            dx1s[ns][0] = (xft[ns + 1][0] - xft[ns][0]) * pdsfilm;
            dx1s[ns][1] = (xft[ns + 1][1] - xft[ns][1]) * pdsfilm;

            rhsT[ns] = rhsT[ns] - (df1s[ns][0] * dx1s[ns][0] + df1s[ns][1] * dx1s[ns][1]);
        }

        for (int ns = 0; ns <= nfilm; ++ns)
        {   //! predicted velocity u* at i
            ufn[ns][0] = (xft[ns][0] - xfn[ns][0]) * pdtfilm;
            ufn[ns][1] = (xft[ns][1] - xfn[ns][1]) * pdtfilm;
        }

        ufn[0][0] = velfix[0] + accfix[0] * dtfilm;    //! velocity of fixed end
        ufn[0][1] = velfix[1] + accfix[1] * dtfilm;

        for (int ns = 0; ns < nfilm; ++ns)
        {   //! (du*1s * du*1s) at i+1/2
            dx1s[ns][0] = (ufn[ns + 1][0] - ufn[ns][0]) * pdsfilm;
            dx1s[ns][1] = (ufn[ns + 1][1] - ufn[ns][1]) * pdsfilm;

            rhsT[ns] -= (dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1]);
        }

        for (int ns = 0; ns < nfilm; ++ns)
        {   //! 0.5*(Dt+)(Dt-)(dx1s*dx1s)^n at i+1/2
            rhsT[ns] += 0.5 * p2dtfilm;

            dx1s[ns][0] = (xfn[ns + 1][0] - xfn[ns][0]) * pdsfilm;
            dx1s[ns][1] = (xfn[ns + 1][1] - xfn[ns][1]) * pdsfilm;

            rhsT[ns] -= (dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1]) * p2dtfilm;

            dx1s[ns][0] = (xfo[ns + 1][0] - xfo[ns][0]) * pdsfilm;
            dx1s[ns][1] = (xfo[ns + 1][1] - xfo[ns][1]) * pdsfilm;

            rhsT[ns] += (dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1]) * 0.5 * p2dtfilm;
        }

        for (int ns = 0; ns < nfilm; ++ns)
        {   //! dx*1s at i+1/2
            dx1s[ns][0] = (xft[ns + 1][0] - xft[ns][0]) * pdsfilm;
            dx1s[ns][1] = (xft[ns + 1][1] - xft[ns][1]) * pdsfilm;
        }

        //! 计算Tension Matrix的三个对角线元素（TAw、TAp和TAe）和右侧项（rhsT）
        //! (dx*1s)*{(d1s)(d1s)(T(dx1*s)]}  at i+1/2
        TAw[0] = 0.0;
        TAp[0] = -(dx1s[0][0] * dx1s[0][0] + dx1s[0][1] * dx1s[0][1]);
        TAe[0] = dx1s[0][0] * dx1s[1][0] + dx1s[0][1] * dx1s[1][1];

        //! application the condition Eq.(25) in Reference JCP Huang WeiXi
        //! and the acceleration boundary condition of fixed end
        rhsT[0] = rhsT[0] - (dx1s[0][0] * (bendf[0][0] + exf[0][0] - accfix[0]) +
            dx1s[0][1] * (bendf[0][1] + exf[0][1] - accfix[1])) * pdsfilm;
        if (isFree == 2)
        {
            rhsT[nfilm-1] = rhsT[nfilm - 1] - (dx1s[nfilm - 1][0] * (bendf[nfilm - 1][0] + exf[nfilm - 1][0] - accfix[0]) +
                dx1s[nfilm - 1][1] * (bendf[nfilm - 1][1] + exf[nfilm - 1][1] - accfix[1])) * pdsfilm;
        }

        for (int ns = 1; ns < nfilm - 1; ++ns)
        {
            TAw[ns] = dx1s[ns][0] * dx1s[ns - 1][0] + dx1s[ns][1] * dx1s[ns - 1][1];
            TAp[ns] = -2.0 * (dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1]);
            TAe[ns] = dx1s[ns][0] * dx1s[ns + 1][0] + dx1s[ns][1] * dx1s[ns + 1][1];
        }

        if (isFree == 1)
        {
            TAw[nfilm - 1] = dx1s[nfilm - 1][0] * dx1s[nfilm - 2][0] + dx1s[nfilm - 1][1] * dx1s[nfilm - 2][1];
            TAp[nfilm - 1] = -3.0 * (dx1s[nfilm - 1][0] * dx1s[nfilm - 1][0] + dx1s[nfilm - 1][1] * dx1s[nfilm - 1][1]);
            TAe[nfilm - 1] = 0.0;
        }
        else if (isFree == 2)
        {
            TAw[nfilm - 1] = dx1s[nfilm - 1][0] * dx1s[nfilm - 2][0] + dx1s[nfilm - 1][1] * dx1s[nfilm - 2][1];
            TAp[nfilm - 1] = -(dx1s[nfilm - 1][0] * dx1s[nfilm - 1][0] + dx1s[nfilm - 1][1] * dx1s[nfilm - 1][1]);
            TAe[nfilm - 1] = 0.0;
        }

        //! 对右侧项 rhsT 中的每个元素乘以 dsfilm * dsfilm，用 trid 函数解张力的 Poisson 方程
        for (int ns = 0; ns < nfilm; ++ns)
        {
            rhsT[ns] = rhsT[ns] * dsfilm * dsfilm;
        }

        vector<double>::iterator First, Endminus;
        First = TAw.begin(), Endminus = TAw.end() - 1;
        vector<double> TAw_temp(First, Endminus);
        First = TAp.begin(), Endminus = TAp.end() - 1;
        vector<double> TAp_temp(First, Endminus);
        First = TAe.begin(), Endminus = TAe.end() - 1;
        vector<double> TAe_temp(First, Endminus);
        vector<double> rhst_temp = rhsT;

        trid(nfilm, TAw_temp, TAp_temp, TAe_temp, tension, rhst_temp);
    }
    else
    {
        cout << "error filament: please input extensible to 1 or 2!!!" << endl;
        exit(1);
    }

    //! Solve Eq.(26) to obtain the filament position at new time step Xn + 1
    //! (Dt + )(Dt - )x = (d1s)[T(dx1s)] + bendf + exf
    for (int k = 0; k < 2; ++k)
    {
        fAww[0][k] = 0.0;
        fAw[0][k] = 0.0;
        fAp[0][k] = 1.0;
        fAe[0][k] = 0.0;
        fAee[0][k] = 0.0;

        //! Calculate boundary condition at fixed end i=0
        bfilm[0][k] = xfn[0][k] + velfix[k] * dtfilm + 0.5 * accfix[k] * dtfilm * dtfilm;

        if (ifixend == 1)
        {
            fAww[1][k] = 0.0;
            fAw[1][k] = -tension[0] * p2dsfilm;    //! -2.0*kbend*p4dsfilm
            fAp[1][k] = (tension[0] + tension[1]) * p2dsfilm + p2dtfilm;    //! +5.0*kbend*p4dsfilm
            fAe[1][k] = -tension[1] * p2dsfilm;    //! -4.0*kbend*p4dsfilm
            fAee[1][k] = 0.0;    //! +kbend*p4dsfilm
            //! Calculate position at i=1
            bfilm[1][k] = exf[1][k] + (2.0 * xfn[1][k] - xfo[1][k]) * p2dtfilm + bendf[1][k];
        }
        else if (ifixend == 2)
        {
            fAww[1][k] = 0.0;
            fAw[1][k] = -tension[0] * p2dsfilm;
            fAp[1][k] = (tension[0] + tension[1]) * p2dsfilm + p2dtfilm;    //! +7.0*kbend*p4dsfilm
            fAe[1][k] = -tension[1] * p2dsfilm;    //! -4.0*kbend*p4dsfilm
            fAee[1][k] = 0.0;    //! +kbend * p4dsfilm
            //! Calculate position at i=1
            bfilm[1][k] = exf[1][k] + (2.0 * xfn[1][k] - xfo[1][k]) * p2dtfilm + bendf[1][k];
        }
        else
        {
            cout << "error filament: please input bounary condition to 1 or 2!!!" << endl;
            exit(1);
        }

        for (int ns = 2; ns <= nfilm - 2; ++ns)
        {
            fAww[ns][k] = 0.0 + kbend * p4dsfilm;
            fAw[ns][k] = -tension[ns - 1] * p2dsfilm - 4.0 * kbend * p4dsfilm;
            fAp[ns][k] = (tension[ns - 1] + tension[ns]) * p2dsfilm + p2dtfilm + 6.0 * kbend * p4dsfilm;
            fAe[ns][k] = -tension[ns] * p2dsfilm - 4.0 * kbend * p4dsfilm;
            fAee[ns][k] = 0.0 + kbend * p4dsfilm;
            bfilm[ns][k] = exf[ns][k] + (2.0 * xfn[ns][k] - xfo[ns][k]) * p2dtfilm;    //! + bendf[ns][k];
        }

        fAww[nfilm - 1][k] = 0.0;    //! + kbend * p4dsfilm;
        fAw[nfilm - 1][k] = -tension[nfilm - 2] * p2dsfilm;    //! - 4.0 * kbend * p4dsfilm;
        fAp[nfilm - 1][k] = (tension[nfilm - 2] + tension[nfilm - 1]) * p2dsfilm + p2dtfilm;    //! + 5.0 * kbend * p4dsfilm;
        fAe[nfilm - 1][k] = -tension[nfilm - 1] * p2dsfilm;    //! - 2.0 * kbend * p4dsfilm;
        fAee[nfilm - 1][k] = 0.0;
        bfilm[nfilm - 1][k] = exf[nfilm - 1][k] + (2.0 * xfn[nfilm - 1][k] - xfo[nfilm - 1][k]) * p2dtfilm + bendf[nfilm - 1][k];

        if (isFree == 1)
        {
            fAww[nfilm][k] = 0.0;    //! + 2.0 * kbend * p4dsfilm;
            fAw[nfilm][k] = -2.0 * tension[nfilm - 1] * p2dsfilm;    //! - 4.0 * kbend * p4dsfilm;
            fAp[nfilm][k] = 2.0 * tension[nfilm - 1] * p2dsfilm + p2dtfilm;    //! + 2.0 * kbend * p4dsfilm;
            fAe[nfilm][k] = 0.0;
            fAee[nfilm][k] = 0.0;

            bfilm[nfilm][k] = exf[nfilm][k] + (2.0 * xfn[nfilm][k] - xfo[nfilm][k]) * p2dtfilm + bendf[nfilm][k];
        }
        else if (isFree == 2)
        {
            fAww[nfilm][k] = 0.0;
            fAw[nfilm][k] = 0.0;
            fAp[nfilm][k] = 1.0;
            fAe[nfilm][k] = 0.0;
            fAee[nfilm][k] = 0.0;

            bfilm[nfilm][k] = xfn[nfilm][k] + velfix[k] * dtfilm + 0.5 * accfix[k] * dtfilm * dtfilm;
        }

        //! ************改为每个列向量**************
        vector<double> fAww_temp(nfilm + 1, 0.0);
        vector<double> fAw_temp(nfilm + 1, 0.0);
        vector<double> fAp_temp(nfilm + 1, 0.0);
        vector<double> fAe_temp(nfilm + 1, 0.0);
        vector<double> fAee_temp(nfilm + 1, 0.0);
        vector<double> xfa_temp(nfilm + 1, 0.0);
        vector<double> bfilm_temp(nfilm + 1, 0.0);
        for (int ntemp = 0; ntemp <= nfilm; ++ntemp)
        {
            fAww_temp[ntemp] = fAww[ntemp][k];
            fAw_temp[ntemp] = fAw[ntemp][k];
            fAp_temp[ntemp] = fAp[ntemp][k];
            fAe_temp[ntemp] = fAe[ntemp][k];
            fAee_temp[ntemp] = fAee[ntemp][k];
            xfa_temp[ntemp] = xfa[ntemp][k];
            bfilm_temp[ntemp] = bfilm[ntemp][k];
        }

        pdma2(nfilm + 1, fAww_temp, fAw_temp, fAp_temp, fAe_temp, fAee_temp, xfa_temp, bfilm_temp);
        for (int ntemp = 0; ntemp <= nfilm; ++ntemp)
        {
            xfa[ntemp][k] = xfa_temp[ntemp];
        }
    }

    eps = 0.0;

    //! Calculate the magnitude of the difference between adjacent nodes
    for (int ns = 0; ns < nfilm; ++ns)
    {
        dx1s[ns][0] = (xfa[ns + 1][0] - xfa[ns][0]) * pdsfilm;
        dx1s[ns][1] = (xfa[ns + 1][1] - xfa[ns][1]) * pdsfilm;

        teps = fabs(dx1s[ns][0] * dx1s[ns][0] + dx1s[ns][1] * dx1s[ns][1] - 1.0);
        if (teps > eps)
        {
            eps = teps;
        }
    }

    //! Update the old and new positions of the nodes
    for (int ns = 0; ns <= nfilm; ++ns)
    {
        xfo[ns][0] = xfn[ns][0];
        xfo[ns][1] = xfn[ns][1];
        xfn[ns][0] = xfa[ns][0];
        xfn[ns][1] = xfa[ns][1];
    }

    //! Calculate the velocity of Lagrangian nodes
    for (int ns = 0; ns <= nfilm; ++ns)
    {
        ufn[ns][0] = (xfn[ns][0] - xfo[ns][0]) * pdtfilm;
        ufn[ns][1] = (xfn[ns][1] - xfo[ns][1]) * pdtfilm;
    }
}

void IBFSI::init_VIV()
{
    double Cylinder_MassRatio = VIV_MassRotio;
    double Cylinder_Diameter_D = Lref;
    double Cylinder_Uniform_U = Uref;

    double Cylinder_Mass;    //! 质量
    double Cylinder_Inertia_Moment;

    double Cylinder_NonDimNatFreq_y;    //! 两个方向无量纲固有频率
    double Cylinder_NonDimNatFreq_Ang;
    double Cylinder_DampRatio_y;    //! 阻尼比
    double Cylinder_DampRatio_Ang;
    Cylinder_NonDimNatFreq_y = VIV_NonDimNatFreq_y;
    Cylinder_NonDimNatFreq_Ang = VIV_NonDimNatFreq_Ang;
    Cylinder_DampRatio_y = VIV_DampRatio_y;
    Cylinder_DampRatio_Ang = VIV_DampRatio_Ang;
    //Cylinder_NonDimNatFreq_y = VIV_NonDimNatFreq_y;
    //Cylinder_NonDimNatFreq_Ang = VIV_NonDimNatFreq_Ang;
    //Cylinder_DampRatio_y = VIV_DampRatio_y;
    //Cylinder_DampRatio_Ang = VIV_DampRatio_Ang;

    vector<double> Cylinder_X(4, 0.0);
    vector<double> Cylinder_Fn(4, 0.0);
    vector<double> Cylinder_Fn_1(4, 0.0);
    vector<double> Cylinder_Fn_2(4, 0.0);
    vector<double> Cylinder_Fn_3(4, 0.0);
    vector<double> Cylinder_Qn(2, 0.0);
    vector<double> Cylinder_Qn_1(2, 0.0);
    vector<double> Cylinder_Qn_2(2, 0.0);
    vector<double> Cylinder_Qn_3(2, 0.0);

    vector<vector<double>> Cylinder_A(4, vector<double>(4, 0.0));
    vector<vector<double>> Cylinder_B(4, vector<double>(2, 0.0));

    Cylinder_Mass = Cylinder_MassRatio * (pi * denIn * Cylinder_Diameter_D * Cylinder_Diameter_D) / 4.0;
    Cylinder_Inertia_Moment = 0.5 * Cylinder_Mass * (Cylinder_Diameter_D / 2) * (Cylinder_Diameter_D / 2);
 
    double Cylinder_K_y = Cylinder_Mass * pow(2 * pi * Cylinder_NonDimNatFreq_y * Cylinder_Uniform_U / Cylinder_Diameter_D, 2);
    double Cylinder_K_Ang = Cylinder_Inertia_Moment * pow(2 * pi * Cylinder_NonDimNatFreq_Ang * Cylinder_Uniform_U / Cylinder_Diameter_D, 2);

    double Cylinder_Damp_y = 2.0 * Cylinder_DampRatio_y * sqrt(Cylinder_K_y * Cylinder_Mass);
    double Cylinder_Damp_Ang = 2.0 * Cylinder_DampRatio_Ang * sqrt(Cylinder_K_Ang * Cylinder_Inertia_Moment);

    Cylinder_A[0][2] = 1.0;
    Cylinder_A[1][3] = 1.0;
    Cylinder_A[2][0] = -Cylinder_K_y / Cylinder_Mass;
    Cylinder_A[3][1] = -Cylinder_K_Ang / Cylinder_Inertia_Moment;
    Cylinder_A[2][2] = -Cylinder_Damp_y / Cylinder_Mass;
    Cylinder_A[3][3] = -Cylinder_Damp_Ang / Cylinder_Inertia_Moment;

    Cylinder_B[2][0] = 1.0 / Cylinder_Mass;
    Cylinder_B[3][1] = 1.0 / Cylinder_Inertia_Moment;

    VIV_X = Cylinder_X;
    VIV_Fn = Cylinder_Fn;
    VIV_Fn_1 = Cylinder_Fn_1;
    VIV_Fn_2 = Cylinder_Fn_2;
    VIV_Fn_3 = Cylinder_Fn_3;
    VIV_Qn = Cylinder_Qn;
    VIV_Qn_1 = Cylinder_Qn_1;
    VIV_Qn_2 = Cylinder_Qn_2;
    VIV_Qn_3 = Cylinder_Qn_3;

    VIV_A = Cylinder_A;
    VIV_B = Cylinder_B;

    VIV_Mass = Cylinder_Mass;    //! 质量
    VIV_Inertia_Moment = Cylinder_Inertia_Moment;

    VIV_K_y = Cylinder_K_y;    //! 有量纲固有频率
    VIV_K_Ang = Cylinder_K_Ang;
    VIV_Damp_y = Cylinder_Damp_y;    //! 有量纲阻尼固有频率
    VIV_Damp_Ang = Cylinder_Damp_Ang;

}

void IBFSI::solve_VIV(double deltat, vector<double> &X, vector<double> &Fn, vector<double> &Fn_1, vector<double> &Fn_2, vector<double> &Fn_3,
    vector<double> &Qn, vector<double> &Qn_1, vector<double> &Qn_2, vector<double> &Qn_3, vector<vector<double>> &A, vector<vector<double>> &B)
{
    vector<double> X_temp(4, 0.0);
    vector<double> F_temp(4, 0.0);
    vector<double> Q_temp(2, 0.0);

    //! t=n
    Fn = add1D_matrix(multiply_two_one(A, X), multiply_two_one(B, Qn));

    for (int i = 0; i < 4; ++i)
    {
        X_temp[i] = X[i] + deltat / 24.0 * (55.0 * Fn[i] - 59.0 * Fn_1[i] + 37.0 * Fn_2[i] - 9.0 * Fn_3[i]);
    }
    Q_temp[0] = 4.0 * Qn[0] - 6.0 * Qn_1[0] + 4.0 * Qn_2[0] - Qn_3[0];
    Q_temp[1] = 4.0 * Qn[1] - 6.0 * Qn_1[1] + 4.0 * Qn_2[1] - Qn_3[1];
    F_temp = add1D_matrix(multiply_two_one(A, X_temp), multiply_two_one(B, Q_temp));

    for (int i = 0; i < 4; ++i)
    {
        X[i] += deltat / 24.0 * (9.0 * F_temp[i] + 19.0 * Fn[i] - 5.0 * Fn_1[i] + Fn_2[i]);
    }

    //! save data 
    Fn_3 = Fn_2;
    Fn_2 = Fn_1;
    Fn_1 = Fn;

    Qn_3 = Qn_2;
    Qn_2 = Qn_1;
    Qn_1 = Qn;
}

//! Subroutine to solve a tridiagonal linear equation system
void trid(int N, vector<double> a, vector<double> b, vector<double> c, vector<double> &x, vector<double> f)
{
    vector<double> temp(N);
    double t;

    c[0] /= b[0];
    f[0] /= b[0];
    for (int k = 1; k < N; ++k)
    {
        t = b[k] - c[k - 1] * a[k];
        c[k] /= t;
        f[k] = (f[k] - f[k - 1] * a[k]) / t;
    }

    temp[N - 1] = f[N - 1];
    for (int k = N - 2; k >= 0; --k)
    {
        temp[k] = f[k] - c[k] * temp[k + 1];
    }

    x = temp;    //! Copy result to output vector x
}

//! 高斯消元法解五对角方程组
void pdma2(int n, vector<double> fAww, vector<double> fAw, vector<double> fAp, vector<double> fAe,
    vector<double> fAee, vector<double> &x, vector<double> frhs)
{
    vector<double> Aww(n), Aw(n), Ap(n), Ae(n), Aee(n), rhs(n);

    Aww = fAww;
    Aw = fAw;
    Ap = fAp;
    Ae = fAe;
    Aee = fAee;
    rhs = frhs;

    //! 消元
    for (int i = 0; i <= n - 3; ++i)
    {
        Ae[i] = Ae[i] / Ap[i];
        Aee[i] = Aee[i] / Ap[i];
        rhs[i] = rhs[i] / Ap[i];

        Ap[i + 1] -= Ae[i] * Aw[i + 1];
        Ae[i + 1] -= Aee[i] * Aw[i + 1];
        rhs[i + 1] -= rhs[i] * Aw[i + 1];

        Aw[i + 2] -= Ae[i] * Aww[i + 2];
        Ap[i + 2] -= Aee[i] * Aww[i + 2];
        rhs[i + 2] -= rhs[i] * Aww[i + 2];
    }

    Ae[n - 2] = Ae[n - 2] / Ap[n - 2];
    Aee[n - 2] = Aee[n - 2] / Ap[n - 2];
    rhs[n - 2] = rhs[n - 2] / Ap[n - 2];

    Ap[n - 1] -= Ae[n - 2] * Aw[n - 1];
    rhs[n - 1] -= rhs[n - 2] * Aw[n - 1];
    rhs[n - 1] = rhs[n - 1] / Ap[n - 1];

    //! 回代
    x[n - 1] = rhs[n - 1];
    x[n - 2] = rhs[n - 2] - Ae[n - 2] * x[n - 1];
    for (int i = n - 3; i >= 0; --i)
    {
        x[i] = rhs[i] - Ae[i] * x[i + 1] - Aee[i] * x[i + 2];
    }
}