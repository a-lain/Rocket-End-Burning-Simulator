#include "../Simulation.hpp"
#include "../Solvers/Gas/ExactRiemannSolver.hpp"
#include "../Solvers/Gas/GasSolvers.hpp"
#include "../Mesh/Mesh.hpp"
#include "../CPGF/Plot2d/LinePlot.hpp"
#include "../Utilities/FormatNumber.hpp"
#include "../Solvers/Gas/ExactSteadySolver.hpp"
#include <iostream>
#include <limits>

using namespace Chemistry;
using namespace Math;
using namespace Utilities;
using namespace CPGF::Plot2d;
using namespace CPGF::AffineSpace;
using namespace CPGF;

void anim_shock_test(const unsigned int i)
{
    double x_0 = 0.5;
    if (i == 4)
    {
        x_0 = 0.3;
    }
    else if (i == 13)
    {
        x_0 = 0.4;
    }
    else if (i == 16)
    {
        x_0 = 0.8;
    }

    Simulation simHLL("test" + std::to_string(i) + "-HLL.sim");
    Simulation simHLLC("test" + std::to_string(i) + "-HLLC.sim");
    Simulation simRoe("test" + std::to_string(i) + "-Roe.sim");
    Simulation simexact("test" + std::to_string(i) + "-exact.sim");

    GasMesh starting_mesh = (*simRoe.gas_mesh)[i];
    Solvers::Gas::ExactRiemannSolver RS(starting_mesh.first_cell->rho,
        starting_mesh.first_cell->v,
        starting_mesh.first_cell->P,
        starting_mesh.last_cell->rho,
        starting_mesh.last_cell->v,
        starting_mesh.last_cell->P,
        starting_mesh.QR.gamma,
        Solvers::Gas::EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE);

    // We determine max and min values.
    double rho_min = std::numeric_limits<double>::max();
    double rho_max = 0;
    double v_min = std::numeric_limits<double>::max();
    double v_max = std::numeric_limits<double>::lowest();
    double P_min = std::numeric_limits<double>::max();
    double P_max = 0;
    for (unsigned int i = 0; i < simRoe.t_array.size(); i++)
    {
        GasMesh mesh = (*simRoe.gas_mesh)[i];
        for (GasMesh::Iterator C = mesh.begin(); C != mesh.end(); ++C)
        {
            if (C->rho < rho_min)
            {
                rho_min = C->rho;
            }
            else if (C->rho > rho_max)
            {
                rho_max = C->rho;
            }
            if (C->v < v_min)
            {
                v_min = C->v;
            }
            else if (C->v > v_max)
            {
                v_max = C->v;
            }
            if (C->P < P_min)
            {
                P_min = C->P;
            }
            else if (C->P > P_max)
            {
                P_max = C->P;
            }
        }
    }

    unsigned int N_points = 300;
    double t = 0;
    double t_video = 5;
    double delta_t = simHLL.t_array[simHLL.t_array.size() - 1] / (t_video*60 - 1);
    for (unsigned int frame = 0; frame < t_video*60; frame++)
    {
        std::vector<double> X;
        std::vector<double> rho_ERS;
        std::vector<double> rho_HLL;
        std::vector<double> rho_HLLC;
        std::vector<double> rho_Roe;
        std::vector<double> rho_exact;
        std::vector<double> v_ERS;
        std::vector<double> v_HLL;
        std::vector<double> v_HLLC;
        std::vector<double> v_Roe;
        std::vector<double> v_exact;
        std::vector<double> P_ERS;
        std::vector<double> P_HLL;
        std::vector<double> P_HLLC;
        std::vector<double> P_Roe;
        std::vector<double> P_exact;
        std::function<double(const double x)> f_rho_HLL = simHLL.rho(t);
        std::function<double(const double x)> f_rho_HLLC = simHLLC.rho(t);
        std::function<double(const double x)> f_rho_Roe = simRoe.rho(t);
        std::function<double(const double x)> f_rho_exact = simexact.rho(t);
        std::function<double(const double x)> f_v_HLL = simHLL.v(t);
        std::function<double(const double x)> f_v_HLLC = simHLLC.v(t);
        std::function<double(const double x)> f_v_Roe = simRoe.v(t);
        std::function<double(const double x)> f_v_exact = simexact.v(t);
        std::function<double(const double x)> f_P_HLL = simHLL.P(t);
        std::function<double(const double x)> f_P_HLLC = simHLLC.P(t);
        std::function<double(const double x)> f_P_Roe = simRoe.P(t);
        std::function<double(const double x)> f_P_exact = simexact.P(t);

        double x = 0;
        double delta_x = 1. / (N_points - 1);

        for (unsigned int i = 0; i < N_points; i++)
        {
            X.push_back(x);
            rho_ERS.push_back(RS.rho(x - x_0, t));
            rho_HLL.push_back(f_rho_HLL(x));
            rho_HLLC.push_back(f_rho_HLLC(x));
            rho_Roe.push_back(f_rho_Roe(x));
            rho_exact.push_back(f_rho_exact(x));
            v_ERS.push_back(RS.v(x - x_0, t));
            v_HLL.push_back(f_v_HLL(x));
            v_HLLC.push_back(f_v_HLLC(x));
            v_Roe.push_back(f_v_Roe(x));
            v_exact.push_back(f_v_exact(x));
            P_ERS.push_back(RS.P(x - x_0, t));
            P_HLL.push_back(f_P_HLL(x));
            P_HLLC.push_back(f_P_HLLC(x));
            P_Roe.push_back(f_P_Roe(x));
            P_exact.push_back(f_P_exact(x));
            x += delta_x;
        }

        LinePlot g_rho_ERS(rho_ERS, X, Color::BLACK, LineWidth::THIN);
        LinePlot g_rho_HLL(rho_HLL, X, Color::BLUE, 1.3*LineWidth::THIN);
        LinePlot g_rho_HLLC(rho_HLLC, X, Color::RED, 1.2*LineWidth::THIN);
        LinePlot g_rho_Roe(rho_Roe, X, Color::GREEN, 1.1*LineWidth::THIN);
        LinePlot g_rho_exact(rho_Roe, X, Color::ORANGE, LineWidth::THIN);
        LinePlot g_v_ERS(v_ERS, X, Color::BLACK, LineWidth::THIN);
        LinePlot g_v_HLL(v_HLL, X, Color::BLUE, 1.3*LineWidth::THIN);
        LinePlot g_v_HLLC(v_HLLC, X, Color::RED, 1.2*LineWidth::THIN);
        LinePlot g_v_Roe(v_Roe, X, Color::GREEN, 1.1*LineWidth::THIN);
        LinePlot g_v_exact(v_Roe, X, Color::ORANGE, LineWidth::THIN);
        LinePlot g_P_ERS(P_ERS, X, Color::BLACK, LineWidth::THIN);
        LinePlot g_P_HLL(P_HLL, X, Color::BLUE, 1.3*LineWidth::THIN);
        LinePlot g_P_HLLC(P_HLLC, X, Color::RED, 1.2*LineWidth::THIN);
        LinePlot g_P_Roe(P_Roe, X, Color::GREEN, 1.1*LineWidth::THIN);
        LinePlot g_P_exact(P_Roe, X, Color::ORANGE, LineWidth::THIN);
        g_rho_ERS.legend = "solución exacta";
        g_rho_HLL.legend = "Godunov-HLL";
        g_rho_HLLC.legend = "Godunov-HLLC";
        g_rho_Roe.legend = "Godunov-Roe";   
        g_rho_exact.legend = "Godunov-exacto";
        g_v_ERS.legend = "solución exacta";
        g_v_HLL.legend = "Godunov-HLL";
        g_v_HLLC.legend = "Godunov-HLLC";
        g_v_Roe.legend = "Godunov-Roe";   
        g_v_exact.legend = "Godunov-exacto";
        g_P_ERS.legend = "solución exacta";
        g_P_HLL.legend = "Godunov-HLL";
        g_P_HLLC.legend = "Godunov-HLLC";
        g_P_Roe.legend = "Godunov-Roe";   
        g_P_exact.legend = "Godunov-exacto";

        Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (u.a.)");
        Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho_g$ (u.a.)");
        Y_axis_rho.set_min_value(rho_min);
        Y_axis_rho.set_max_value(rho_max);
        Graphic G_rho;
        G_rho.show_legend = true;
        G_rho.legend_position = LegendPosition::ABOVE;
        G_rho.add(&g_rho_ERS, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_HLL, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_HLLC, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_Roe, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_exact, &Y_axis_rho, &X_axis_rho);
        Scene2d S_rho = G_rho.render_to_scene();
        Text text(Point2d(0, Axis::X_MAX*0.6 + Graphic::LEGEND_MARGIN + 6*Graphic::LEGEND_VERTICAL_DISPLACEMENT_PER_LINE),
            "$t = " + Utilities::format_number(t, true, 3, false, -10, 10) + "$ s", Color::BLACK, TextAlignment::LEFT);
        S_rho += text;
        char buffer[10];
        sprintf(buffer, "%04d", frame);
        S_rho.render("/home/andres/Animaciones/" + std::to_string(i) + "/rho/" + std::string(buffer) + ".png", 250);

        Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (u.a.)");
        Axis Y_axis_v(AxisType::VERTICAL, "$v$ (u.a.)");
        Y_axis_v.set_min_value(v_min);
        Y_axis_v.set_max_value(v_max);
        Graphic G_v;
        G_v.show_legend = true;
        G_v.legend_position = LegendPosition::ABOVE;
        G_v.add(&g_v_ERS, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_HLL, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_HLLC, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_Roe, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_exact, &Y_axis_v, &X_axis_v);
        Scene2d S_v = G_v.render_to_scene();
        text = Text(Point2d(0, Axis::X_MAX*0.6 + Graphic::LEGEND_MARGIN + 6*Graphic::LEGEND_VERTICAL_DISPLACEMENT_PER_LINE),
            "$t = " + Utilities::format_number(t, true, 3, false, -10, 10) + "$ s", Color::BLACK, TextAlignment::LEFT);
        S_v += text;
        sprintf(buffer, "%04d", frame);
        S_v.render("/home/andres/Animaciones/" + std::to_string(i) + "/v/" + std::string(buffer) + ".png", 250);

        Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (u.a.)");
        Axis Y_axis_P(AxisType::VERTICAL, "$P$ (u.a.)");
        Y_axis_P.set_min_value(P_min);
        Y_axis_P.set_max_value(P_max);
        Graphic G_P;
        G_P.show_legend = true;
        G_P.legend_position = LegendPosition::ABOVE;
        G_P.add(&g_P_ERS, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_HLL, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_HLLC, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_Roe, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_exact, &Y_axis_P, &X_axis_P);
        Scene2d S_P = G_P.render_to_scene();
        text = Text(Point2d(0, Axis::X_MAX*0.6 + Graphic::LEGEND_MARGIN + 6*Graphic::LEGEND_VERTICAL_DISPLACEMENT_PER_LINE),
            "$t = " + Utilities::format_number(t, true, 3, false, -10, 10) + "$ s", Color::BLACK, TextAlignment::LEFT);
        S_P += text;
        sprintf(buffer, "%04d", frame);
        S_P.render("/home/andres/Animaciones/" + std::to_string(i) + "/P/" + std::string(buffer) + ".png", 250);


        std::cout << frame << " / " << t_video*60 << "\n";

        t += delta_t;
    }
}

void anim19()
{
    Simulation sim("testB1-Roe.sim");

    // // We fix QR.
    // for (unsigned int i = 0; i < sim.t_array.size(); i++)
    // {
    //     GasMesh gas_mesh = (*sim.gas_mesh)[i];
    //     SolidMesh solid_mesh = (*sim.solid_mesh)[i];
    //     for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
    //     {
    //         C->QR = &sim.QR;
    //     }
    //     for (SolidMesh::Iterator C = solid_mesh.begin(); C != solid_mesh.end(); ++C)
    //     {
    //         C->QR = &sim.QR;
    //     }
    //     (*sim.gas_mesh)[i] = gas_mesh;
    //     (*sim.solid_mesh)[i] = solid_mesh;
    // }

    unsigned int N_points = 300;
    double t = 0;
    double t_video = sim.t_array[sim.t_array.size() - 1];
    double delta_t = sim.t_array[sim.t_array.size() - 1] / (t_video*60 - 1);
    for (unsigned int frame = 0; frame < t_video*60; frame++)
    {
        std::vector<double> X;
        std::vector<double> X_gas;
        std::vector<double> rho_ERS;
        std::vector<double> rho_Roe;
        std::vector<double> v_ERS;
        std::vector<double> v_Roe;
        std::vector<double> P_ERS;
        std::vector<double> P_Roe;
        std::vector<double> T_ERS;
        std::vector<double> T_Roe;
        std::function<double(const double x)> f_rho_Roe = sim.rho(t);
        std::function<double(const double x)> f_v_Roe = sim.v(t);
        std::function<double(const double x)> f_P_Roe = sim.P(t);
        std::function<double(const double x)> f_T_Roe = sim.T(t);

        double x = 0;
        double delta_x = 0.3 / (N_points - 1);

        for (unsigned int i = 0; i < N_points; i++)
        {
            X.push_back(x);
            T_Roe.push_back(f_T_Roe(x));
            if (x > sim.x_q(t))
            {
                X_gas.push_back(x);
                T_ERS.push_back(1735.29);
                rho_ERS.push_back(2012.6);
                rho_Roe.push_back(f_rho_Roe(x));
                v_ERS.push_back(0);
                v_Roe.push_back(f_v_Roe(x));
                P_ERS.push_back(558.79e6);
                P_Roe.push_back(f_P_Roe(x));
            }
            x += delta_x;
        }

        LinePlot g_rho_ERS(rho_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_rho_Roe(rho_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_v_ERS(v_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_v_Roe(v_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_P_ERS(P_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_P_Roe(P_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_T_ERS(T_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_T_Roe(T_Roe, X, Color::GREEN, LineWidth::THIN);
        g_rho_ERS.legend = "aproximación estacionaria";
        g_rho_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";
        g_v_ERS.legend = "aproximación estacionaria";
        g_v_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_P_ERS.legend = "aproximación estacionaria";
        g_P_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_T_ERS.legend = "aproximación estacionaria";
        g_T_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s"; 

        Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
        Y_axis_rho.set_min_value(0);
        Y_axis_rho.set_max_value(2090);
        X_axis_rho.set_min_value(0);
        Graphic G_rho;
        G_rho.show_legend = true;
        G_rho.legend_position = LegendPosition::ABOVE;
        G_rho.add(&g_rho_ERS, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_Roe, &Y_axis_rho, &X_axis_rho);
        Scene2d S_rho = G_rho.render_to_scene();
        char buffer[10];
        sprintf(buffer, "%04d", frame);
        S_rho.render("/home/andres/Animaciones/B1/rho/" + std::string(buffer) + ".png", 250);

        Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
        Y_axis_v.set_min_value(-5);
        Y_axis_v.set_max_value(5);
        X_axis_v.set_min_value(0);
        Graphic G_v;
        G_v.show_legend = true;
        G_v.legend_position = LegendPosition::ABOVE;
        G_v.add(&g_v_ERS, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_Roe, &Y_axis_v, &X_axis_v);
        Scene2d S_v = G_v.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_v.render("/home/andres/Animaciones/B1/v/" + std::string(buffer) + ".png", 250);

        Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
        Y_axis_P.set_min_value(0);
        Y_axis_P.set_max_value(5.6e8);
        X_axis_P.set_min_value(0);
        Graphic G_P;
        G_P.show_legend = true;
        G_P.legend_position = LegendPosition::ABOVE;
        G_P.add(&g_P_ERS, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_Roe, &Y_axis_P, &X_axis_P);
        Scene2d S_P = G_P.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_P.render("/home/andres/Animaciones/B1/P/" + std::string(buffer) + ".png", 250);

        Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
        Y_axis_T.set_min_value(0);
        Y_axis_T.set_max_value(2000);
        Graphic G_T;
        G_T.show_legend = true;
        G_T.legend_position = LegendPosition::ABOVE;
        G_T.add(&g_T_ERS, &Y_axis_T, &X_axis_T);
        G_T.add(&g_T_Roe, &Y_axis_T, &X_axis_T);
        Scene2d S_T = G_T.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_T.render("/home/andres/Animaciones/B1/T/" + std::string(buffer) + ".png", 250);

        std::cout << frame << " / " << t_video*60 << "\n";

        t += delta_t;
    }
}

void anim20(const unsigned int i_start, const unsigned int i_end)
{
    Simulation sim("testB6-Roe.sim");
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    // We fix QR.
    for (unsigned int i = 0; i < sim.t_array.size(); i++)
    {
        GasMesh gas_mesh = (*sim.gas_mesh)[i];
        SolidMesh solid_mesh = (*sim.solid_mesh)[i];
        gas_mesh.QR = QR;
        solid_mesh.QR = QR;
        (*sim.gas_mesh)[i] = gas_mesh;
        (*sim.solid_mesh)[i] = solid_mesh;
    }

    unsigned int N_points = 300;
    double t_video = sim.t_array[sim.t_array.size() - 1];
    double delta_t = sim.t_array[sim.t_array.size() - 1] / (t_video*60 - 1);
    double t = i_start * delta_t;
    for (unsigned int frame = i_start; frame < i_end; frame++)
    {
        std::vector<double> X;
        std::vector<double> X_gas;
        std::vector<double> rho_ERS;
        std::vector<double> rho_Roe;
        std::vector<double> v_ERS;
        std::vector<double> v_Roe;
        std::vector<double> P_ERS;
        std::vector<double> P_Roe;
        std::vector<double> T_ERS;
        std::vector<double> T_Roe;
        std::function<double(const double x)> f_rho_Roe = sim.rho(t);
        std::function<double(const double x)> f_v_Roe = sim.v(t);
        std::function<double(const double x)> f_P_Roe = sim.P(t);
        std::function<double(const double x)> f_T_Roe = sim.T(t);

        double x = 0;
        double delta_x = 0.3 / (N_points - 1);

        for (unsigned int i = 0; i < N_points; i++)
        {
            X.push_back(x);
            T_Roe.push_back(f_T_Roe(x));
            if (x > sim.x_q(t))
            {
                X_gas.push_back(x);
                T_ERS.push_back(1423.68);
                P_ERS.push_back(101325);
                rho_ERS.push_back(0.444820);
                v_ERS.push_back(37.892452);
                rho_Roe.push_back(f_rho_Roe(x));
                v_Roe.push_back(f_v_Roe(x));
                P_Roe.push_back(f_P_Roe(x));
            }
            x += delta_x;
        }

        LinePlot g_rho_ERS(rho_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_rho_Roe(rho_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_v_ERS(v_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_v_Roe(v_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_P_ERS(P_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_P_Roe(P_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_T_ERS(T_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_T_Roe(T_Roe, X, Color::GREEN, LineWidth::THIN);
        g_rho_ERS.legend = "aproximación estacionaria";
        g_rho_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";
        g_v_ERS.legend = "aproximación estacionaria";
        g_v_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_P_ERS.legend = "aproximación estacionaria";
        g_P_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_T_ERS.legend = "aproximación estacionaria";
        g_T_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s"; 

        Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
        Y_axis_rho.set_min_value(0);
        Y_axis_rho.set_max_value(0.6);
        X_axis_rho.set_min_value(0);
        Graphic G_rho;
        G_rho.show_legend = true;
        G_rho.legend_position = LegendPosition::ABOVE;
        G_rho.add(&g_rho_ERS, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_Roe, &Y_axis_rho, &X_axis_rho);
        Scene2d S_rho = G_rho.render_to_scene();
        char buffer[10];
        sprintf(buffer, "%04d", frame);
        S_rho.render("/home/andres/Animaciones/B6/rho/" + std::string(buffer) + ".png", 250);

        Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
        Y_axis_v.set_min_value(0);
        Y_axis_v.set_max_value(43);
        X_axis_v.set_min_value(0);
        Graphic G_v;
        G_v.show_legend = true;
        G_v.legend_position = LegendPosition::ABOVE;
        G_v.add(&g_v_ERS, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_Roe, &Y_axis_v, &X_axis_v);
        Scene2d S_v = G_v.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_v.render("/home/andres/Animaciones/B6/v/" + std::string(buffer) + ".png", 250);

        Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
        Y_axis_P.set_min_value(0);
        Y_axis_P.set_max_value(103000);
        X_axis_P.set_min_value(0);
        Graphic G_P;
        G_P.show_legend = true;
        G_P.legend_position = LegendPosition::ABOVE;
        G_P.add(&g_P_ERS, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_Roe, &Y_axis_P, &X_axis_P);
        Scene2d S_P = G_P.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_P.render("/home/andres/Animaciones/B6/P/" + std::string(buffer) + ".png", 250);

        Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
        Y_axis_T.set_min_value(290);
        Y_axis_T.set_max_value(1600);
        X_axis_T.set_min_value(0);
        Graphic G_T;
        G_T.show_legend = true;
        G_T.legend_position = LegendPosition::ABOVE;
        G_T.add(&g_T_ERS, &Y_axis_T, &X_axis_T);
        G_T.add(&g_T_Roe, &Y_axis_T, &X_axis_T);
        Scene2d S_T = G_T.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_T.render("/home/andres/Animaciones/B6/T/" + std::string(buffer) + ".png", 250);

        std::cout << frame << " / " << t_video*60 << "\n";

        t += delta_t;
    }
}

void anim21(const unsigned int i_start, const unsigned int i_end)
{
    Simulation sim("testB11-Roe.sim");

    // // We fix QR.
    // for (unsigned int i = 0; i < sim.t_array.size(); i++)
    // {
    //     GasMesh gas_mesh = (*sim.gas_mesh)[i];
    //     SolidMesh solid_mesh = (*sim.solid_mesh)[i];
    //     for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
    //     {
    //         C->QR = &sim.QR;
    //     }
    //     for (SolidMesh::Iterator C = solid_mesh.begin(); C != solid_mesh.end(); ++C)
    //     {
    //         C->QR = &sim.QR;
    //     }
    //     (*sim.gas_mesh)[i] = gas_mesh;
    //     (*sim.solid_mesh)[i] = solid_mesh;
    // }

    double L = 0.35;
    double A_c = M_PI*0.05*0.05;
    double A_e = A_c / 10;
    double x_c = 0.3;
    double x_e = 0.33;
    double L_t = x_e - x_c;
    double x_m = (x_c + x_e) / 2;
    auto A = [=](double x)
    {
        if (x <= x_c)
        {
            return A_c;
        }
        else if (x > x_c && x < x_m)
        {
            return -2*(A_c - A_e)/L_t/L_t*(x - x_c)*(x - x_c) + A_c;
        }
        else if (x > x_m && x < x_e)
        {
            return 2*(A_c - A_e)*((x-x_c)/L_t - 1)*((x-x_c)/L_t - 1) + A_e;
        }
        else
        {
            return A_e;
        }
    };
    Solvers::Gas::ExactSteadySolver SS(0.63813475, 145411, 1424.18, 0.056426, A_c, 1.218878, Solvers::Gas::SolutionType::SUBSONIC);

    unsigned int N_points = 300;
    double t_video = sim.t_array[sim.t_array.size() - 1];
    double delta_t = sim.t_array[sim.t_array.size() - 1] / (t_video*60 - 1);
    double t = i_start * delta_t;
    for (unsigned int frame = i_start; frame < i_end; frame++)
    {
        std::vector<double> X;
        std::vector<double> X_gas;
        std::vector<double> rho_ERS;
        std::vector<double> rho_Roe;
        std::vector<double> v_ERS;
        std::vector<double> v_Roe;
        std::vector<double> P_ERS;
        std::vector<double> P_Roe;
        std::vector<double> T_ERS;
        std::vector<double> T_Roe;
        std::function<double(const double x)> f_rho_Roe = sim.rho(t);
        std::function<double(const double x)> f_v_Roe = sim.v(t);
        std::function<double(const double x)> f_P_Roe = sim.P(t);
        std::function<double(const double x)> f_T_Roe = sim.T(t);

        double x = 0;
        double delta_x = 0.35 / (N_points - 1);

        for (unsigned int i = 0; i < N_points; i++)
        {
            X.push_back(x);
            T_Roe.push_back(f_T_Roe(x));
            if (x > sim.x_q(t))
            {
                X_gas.push_back(x);
                T_ERS.push_back(SS.T(A(x)));
                P_ERS.push_back(SS.P(A(x)));
                rho_ERS.push_back(SS.rho(A(x)));
                v_ERS.push_back(SS.v(A(x)));
                rho_Roe.push_back(f_rho_Roe(x));
                v_Roe.push_back(f_v_Roe(x));
                P_Roe.push_back(f_P_Roe(x));
            }
            x += delta_x;
        }

        LinePlot g_rho_ERS(rho_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_rho_Roe(rho_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_v_ERS(v_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_v_Roe(v_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_P_ERS(P_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_P_Roe(P_Roe, X_gas, Color::GREEN, LineWidth::THIN);
        LinePlot g_T_ERS(T_ERS, X_gas, Color::BLACK, LineWidth::THIN);
        LinePlot g_T_Roe(T_Roe, X, Color::GREEN, LineWidth::THIN);
        g_rho_ERS.legend = "aproximación estacionaria";
        g_rho_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";
        g_v_ERS.legend = "aproximación estacionaria";
        g_v_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_P_ERS.legend = "aproximación estacionaria";
        g_P_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s";  
        g_T_ERS.legend = "aproximación estacionaria";
        g_T_Roe.legend = "$t = " + Utilities::format_number(t) + "$ s"; 

        Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
        Y_axis_rho.set_min_value(0.4);
        Y_axis_rho.set_max_value(0.8);
        X_axis_rho.set_min_value(0);
        Graphic G_rho;
        G_rho.show_legend = true;
        G_rho.legend_position = LegendPosition::ABOVE;
        G_rho.add(&g_rho_ERS, &Y_axis_rho, &X_axis_rho);
        G_rho.add(&g_rho_Roe, &Y_axis_rho, &X_axis_rho);
        Scene2d S_rho = G_rho.render_to_scene();
        char buffer[10];
        sprintf(buffer, "%04d", frame);
        S_rho.render("/home/andres/Animaciones/B11/rho/" + std::string(buffer) + ".png", 250);

        Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
        Y_axis_v.set_min_value(0);
        Y_axis_v.set_max_value(430);
        X_axis_v.set_min_value(0);
        Graphic G_v;
        G_v.show_legend = true;
        G_v.legend_position = LegendPosition::ABOVE;
        G_v.add(&g_v_ERS, &Y_axis_v, &X_axis_v);
        G_v.add(&g_v_Roe, &Y_axis_v, &X_axis_v);
        Scene2d S_v = G_v.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_v.render("/home/andres/Animaciones/B11/v/" + std::string(buffer) + ".png", 250);

        Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
        Y_axis_P.set_min_value(1e5);
        Y_axis_P.set_max_value(1.6e5);
        X_axis_P.set_min_value(0);
        Graphic G_P;
        G_P.show_legend = true;
        G_P.legend_position = LegendPosition::ABOVE;
        G_P.add(&g_P_ERS, &Y_axis_P, &X_axis_P);
        G_P.add(&g_P_Roe, &Y_axis_P, &X_axis_P);
        Scene2d S_P = G_P.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_P.render("/home/andres/Animaciones/B11/P/" + std::string(buffer) + ".png", 250);

        Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
        Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
        Y_axis_T.set_min_value(290);
        Y_axis_T.set_max_value(1500);
        X_axis_T.set_min_value(0);
        Graphic G_T;
        G_T.show_legend = true;
        G_T.legend_position = LegendPosition::ABOVE;
        G_T.add(&g_T_ERS, &Y_axis_T, &X_axis_T);
        G_T.add(&g_T_Roe, &Y_axis_T, &X_axis_T);
        Scene2d S_T = G_T.render_to_scene();
        sprintf(buffer, "%04d", frame);
        S_T.render("/home/andres/Animaciones/B11/T/" + std::string(buffer) + ".png", 250);

        std::cout << frame << " / " << t_video*60 << "\n";

        t += delta_t;
    }
}
