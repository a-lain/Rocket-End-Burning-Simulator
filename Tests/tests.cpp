#include "../Simulation.hpp"
#include "../Solvers/Gas/GasSolvers.hpp"
#include "../Solvers/Gas/ExactRiemannSolver.hpp"
#include "../Solvers/Gas/ExactSteadySolver.hpp"
#include "../Solvers/Solid/SolidSolvers.hpp"
#include "../Chemistry/Reaction.hpp"
#include "../Math/Vector.hpp"
#include "../CPGF/Plot2d/LinePlot.hpp"
#include <iostream>
#include <cmath>

using namespace Chemistry;
using namespace Math;
using namespace Utilities;
using namespace CPGF::Plot2d;
using namespace CPGF;

void test0()
{
    GasBoundaryConditions BC(GasBoundaryConditionsType::WALL, GasBoundaryConditionsType::WALL);
    double rho_L = 1;
    double v_L = 0;
    double P_L = 1;
    double rho_R = 1;
    double v_R = 0;
    double P_R = 1;
    double t_end = 1000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;
    
    SolidGasReaction QR(0, 0, 0, 2.5, 1, 0, 101325, 0, 0);

    Simulation sim1("test0-HLL", QR,
                0, 1, N_cells, false, 1, 0.9, 0, BC,
                [=](double x){return (x < 0.5) ? rho_L : rho_R;},
                [=](double x){return (x < 0.5) ? v_L : v_R;},
                [=](double x){return (x < 0.5) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLL);
    sim1.simulate_until(t_end);
    
    Simulation sim2("test0-HLLC", QR,
                0, 1, N_cells, false, 1, 0.9, 0, BC,
                [=](double x){return (x < 0.5) ? rho_L : rho_R;},
                [=](double x){return (x < 0.5) ? v_L : v_R;},
                [=](double x){return (x < 0.5) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLLC);
    sim2.simulate_until(t_end);
    
    Simulation sim3("test0-Roe", QR,
                0, 1, N_cells, false, 1, 0.9, 0, BC,
                [=](double x){return (x < 0.5) ? rho_L : rho_R;},
                [=](double x){return (x < 0.5) ? v_L : v_R;},
                [=](double x){return (x < 0.5) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    
    Simulation sim4("test0-exact", QR,
                0, 1, N_cells, false, 1, 0.9, 0, BC,
                [=](double x){return (x < 0.5) ? rho_L : rho_R;},
                [=](double x){return (x < 0.5) ? v_L : v_R;},
                [=](double x){return (x < 0.5) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::exact);
    sim4.simulate_until(t_end);

    // The exact solution.
    std::vector<double> X;
    std::vector<double> rho;
    std::vector<double> v;
    std::vector<double> P;
    double delta_x = 1. / (N_points - 1);
    double x = 0;
    Solvers::Gas::ExactRiemannSolver RS(rho_L, v_L, P_L, rho_R, v_R, P_R, QR.gamma,
        Solvers::Gas::EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE);
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        rho.push_back(RS.rho(x - 0.5, t_end));
        v.push_back(RS.v(x - 0.5, t_end));
        P.push_back(RS.P(x - 0.5, t_end));
        x += delta_x;
    }
    LinePlot g10(rho, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(v, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    
    // The approximated solutions.
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g21(mesh.v(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    AveragePlot g22(mesh.v(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g23(mesh.v(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    mesh = sim4.instant_gas_mesh;
    AveragePlot g14(mesh.rho(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g24(mesh.v(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g34(mesh.P(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    g14.legend = "Godunov-exacto";
    g24.legend = "Godunov-exacto";
    g34.legend = "Godunov-exacto";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g$ (u.a.)");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.add(&g14, &Y_axis1, &X_axis1);
    Y_axis1.set_min_value(0);
    G1.render_to_scene().render("T1.pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T1small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis2(AxisType::VERTICAL, "$v$ (u.a.)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.add(&g24, &Y_axis2, &X_axis2);
    Y_axis2.set_min_value(-1);
    Y_axis2.set_max_value(1);
    G2.render_to_scene().render("T2.pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T2small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (u.a.)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.add(&g34, &Y_axis3, &X_axis3);
    Y_axis3.set_min_value(0);
    G3.render_to_scene().render("T3.pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T3small.pdf");
}

void shock_test(
    const double rho_L, const double v_L, const double P_L,
    const double rho_R, const double v_R, const double P_R,
    const unsigned int i, const double x_0,
    const double t_end, unsigned int N_cells, unsigned int N_points)
{
    GasBoundaryConditions BC(GasBoundaryConditionsType::FREE, GasBoundaryConditionsType::FREE);
    
    SolidGasReaction QR(0, 0, 0, 2.5, 1, 0, 101325, 0, 0);

    Simulation sim1("test" + std::to_string(i) + "-HLL", QR,
                0, 1, N_cells, false, 1, 0.9, 300, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLL);
    sim1.simulate_until(t_end);
    sim1.write_to_file("test" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("test" + std::to_string(i) + "-HLLC", QR,
                0, 1, N_cells, false, 1, 0.9, 300, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLLC);
    sim2.simulate_until(t_end);
    sim2.write_to_file("test" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, 1, N_cells, false, 1, 0.9, 300, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");
    
    Simulation sim4("test" + std::to_string(i) + "-exact", QR,
                0, 1, N_cells, false, 1, 0.9, 300, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::exact);
    sim4.simulate_until(t_end);
    sim4.write_to_file("test" + std::to_string(i) + "-exact.sim");

    // The exact solution.
    std::vector<double> X;
    std::vector<double> rho;
    std::vector<double> v;
    std::vector<double> P;
    double delta_x = 1. / (N_points - 1);
    double x = 0;
    Solvers::Gas::ExactRiemannSolver RS(rho_L, v_L, P_L, rho_R, v_R, P_R, QR.gamma,
        Solvers::Gas::EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE);
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        rho.push_back(RS.rho(x - x_0, t_end));
        v.push_back(RS.v(x - x_0, t_end));
        P.push_back(RS.P(x - x_0, t_end));
        x += delta_x;
    }
    LinePlot g10(rho, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(v, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    
    // The approximated solutions.
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g21(mesh.v(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    AveragePlot g22(mesh.v(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g23(mesh.v(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    mesh = sim4.instant_gas_mesh;
    AveragePlot g14(mesh.rho(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g24(mesh.v(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g34(mesh.P(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    g14.legend = "Godunov-exacto";
    g24.legend = "Godunov-exacto";
    g34.legend = "Godunov-exacto";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g$ (u.a.)");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.add(&g14, &Y_axis1, &X_axis1);
    if (i != 19 && i != 22)
    {
        Y_axis1.set_min_value(0);
    }
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis2(AxisType::VERTICAL, "$v$ (u.a.)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.add(&g24, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (u.a.)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.add(&g34, &Y_axis3, &X_axis3);
    Y_axis3.set_min_value(0);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");
}

void shock_test_adaptative(
    const double rho_L, const double v_L, const double P_L,
    const double rho_R, const double v_R, const double P_R,
    const unsigned int i, const double x_0,
    const double t_end, unsigned int N_cells, unsigned int N_points, unsigned int adaptive_refinement_period)
{
    GasBoundaryConditions BC(GasBoundaryConditionsType::FREE, GasBoundaryConditionsType::FREE);
    
    SolidGasReaction QR(0, 0, 0, 2.5, 1, 0, 101325, 0, 0);

    Simulation sim1("testa" + std::to_string(i) + "-HLL", QR,
                0, 1, N_cells, true, 1, 0.9, 100, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLL,
                [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
                [] (const GasCell& C, const double t){return 0.;},
                1, 0.01, 0.001, 1./50, 1./10000, 1./2000);
    sim1.adaptive_refinement_period = adaptive_refinement_period;
    sim1.simulate_until(t_end);
    sim1.write_to_file("testa" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("testa" + std::to_string(i) + "-HLLC", QR,
                0, 1, N_cells, true, 1, 0.9, 100, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::HLLC,
                [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
                [] (const GasCell& C, const double t){return 0.;},
                1, 0.01, 0.001, 1./50, 1./10000, 1./2000);
    sim2.adaptive_refinement_period = adaptive_refinement_period;
    sim2.simulate_until(t_end); 
    sim2.write_to_file("testa" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, 1, N_cells, true, 1, 0.9, 100, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::Roe,
                [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
                [] (const GasCell& C, const double t){return 0.;},
                1, 0.01, 0.001, 1./50, 1./10000, 1./2000);
    sim3.adaptive_refinement_period = adaptive_refinement_period;
    sim3.simulate_until(t_end);
    sim3.write_to_file("testa" + std::to_string(i) + "-Roe.sim");
    
    Simulation sim4("testa" + std::to_string(i) + "-exact", QR,
                0, 1, N_cells, true, 1, 0.9, 100, BC,
                [=](double x){return (x < x_0) ? rho_L : rho_R;},
                [=](double x){return (x < x_0) ? v_L : v_R;},
                [=](double x){return (x < x_0) ? P_L : P_R;},
                [](double x){return 1.;},
                Solvers::Gas::exact,
                [] (std::function<double(double x)> f, const double a, const double b)
                {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
                [] (const GasCell& C, const double t){return 0.;},
                1, 0.01, 0.001, 1./50, 1./10000, 1./2000);
    sim4.adaptive_refinement_period = adaptive_refinement_period;
    sim4.simulate_until(t_end);
    sim4.write_to_file("testa" + std::to_string(i) + "-exact.sim");

    // The exact solution.
    std::vector<double> X;
    std::vector<double> rho;
    std::vector<double> v;
    std::vector<double> P;
    double delta_x = 1. / (N_points - 1);
    double x = 0;
    Solvers::Gas::ExactRiemannSolver RS(rho_L, v_L, P_L, rho_R, v_R, P_R, QR.gamma,
        Solvers::Gas::EXACT_RIEMANN_SOLVER_RELATIVE_TOLERANCE);
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        rho.push_back(RS.rho(x - x_0, t_end));
        v.push_back(RS.v(x - x_0, t_end));
        P.push_back(RS.P(x - x_0, t_end));
        x += delta_x;
    }
    LinePlot g10(rho, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(v, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    
    // The approximated solutions.
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g21(mesh.v(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.3*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g22(mesh.v(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.2*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g23(mesh.v(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, 1.1*LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    mesh = sim4.instant_gas_mesh;
    AveragePlot g14(mesh.rho(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g24(mesh.v(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    AveragePlot g34(mesh.P(), mesh.x_partition(), Color::ORANGE, LineWidth::THIN);
    g14.legend = "Godunov-exacto";
    g24.legend = "Godunov-exacto";
    g34.legend = "Godunov-exacto";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g$ (u.a.)");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.add(&g14, &Y_axis1, &X_axis1);
    if (i != 19 && i != 22)
    {
        Y_axis1.set_min_value(0);
    }
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis2(AxisType::VERTICAL, "$v$ (u.a.)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.add(&g24, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (u.a.)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (u.a.)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.add(&g34, &Y_axis3, &X_axis3);
    Y_axis3.set_min_value(0);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Graphic* GMP1 = sim1.mesh_plot();
    GMP1->render_to_scene().render("M" + std::to_string(i) + "P1.pdf");
    GMP1->axes[0]->scale = 0.5;
    GMP1->axes[0]->N_major_ticks = 5;
    GMP1->axes[1]->scale = 0.5;
    GMP1->axes[1]->N_major_ticks = 5;
    GMP1->render_to_scene().render("M" + std::to_string(i) + "P1small.pdf");
    Graphic* GMP2 = sim2.mesh_plot();
    GMP2->render_to_scene().render("M" + std::to_string(i) + "P2.pdf");
    GMP2->axes[0]->scale = 0.5;
    GMP2->axes[0]->N_major_ticks = 5;
    GMP2->axes[1]->scale = 0.5;
    GMP2->axes[1]->N_major_ticks = 5;
    GMP2->render_to_scene().render("M" + std::to_string(i) + "P2small.pdf");
    Graphic* GMP3 = sim3.mesh_plot();
    GMP3->render_to_scene().render("M" + std::to_string(i) + "P3.pdf");
    GMP3->axes[0]->scale = 0.5;
    GMP3->axes[0]->N_major_ticks = 5;
    GMP3->axes[1]->scale = 0.5;
    GMP3->axes[1]->N_major_ticks = 5;
    GMP3->render_to_scene().render("M" + std::to_string(i) + "P3small.pdf");
    Graphic* GMP4 = sim4.mesh_plot();
    GMP4->render_to_scene().render("M" + std::to_string(i) + "P4.pdf");
    GMP4->axes[0]->scale = 0.5;
    GMP4->axes[0]->N_major_ticks = 5;
    GMP4->axes[1]->scale = 0.5;
    GMP4->axes[1]->N_major_ticks = 5;
    GMP4->render_to_scene().render("M" + std::to_string(i) + "P4small.pdf");
}

void test1()
{
    shock_test(1, 0.75, 1, 0.125, 0, 0.1,
        4, 0.3, 0.2, 200, 300);
}

void test1a()
{
    shock_test_adaptative(1, 0.75, 1, 0.125, 0, 0.1,
        4, 0.3, 0.2, 10000, 300, 1);
}

void test2()
{
    shock_test(1, -2, 3, 1, 2, 3, 
        7, 0.5, 0.15, 200, 300);
}

void test2a()
{
    shock_test_adaptative(1, -2, 3, 1, 2, 3, 
        7, 0.5, 0.15, 10000, 300, 1);
}

void test3()
{
    shock_test(1, 0, 1000, 1, 0, 0.01, 
        10, 0.5, 0.011, 200, 300);
}

void test3a()
{
    shock_test_adaptative(1, 0, 1000, 1, 0, 0.01, 
        10, 0.5, 0.011, 10000, 300, 1);
}

void test4()
{
    shock_test(5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950, 
        13, 0.4, 0.035, 200, 300);
}

void test4a()
{
    shock_test_adaptative(5.99924, 19.5975, 460.894, 5.99242, -6.19633, 46.0950, 
        13, 0.4, 0.035, 10000, 300, 1);
}

void test5()
{
    shock_test(1, -19.59745, 1000, 1, -19.59745, 0.01, 
        16, 0.8, 0.012, 200, 300);
}

void test5a()
{
    shock_test_adaptative(1, -19.59745, 1000, 1, -19.59745, 0.01, 
        16, 0.8, 0.012, 10000, 300, 1);
}

void test6()
{
    shock_test(1.4, 0, 1, 1, 0, 1, 
        19, 0.5, 2, 200, 300);
}

void test6a()
{
    shock_test_adaptative(1.4, 0, 1, 1, 0, 1, 
        19, 0.5, 2, 10000, 300, 1);
}

void test7()
{
    shock_test(1.4, 0.1, 1, 1, 0.1, 1, 
        22, 0.5, 2, 200, 300);
}

void test7a()
{
    shock_test_adaptative(1.4, 0.1, 1, 1, 0.1, 1, 
        22, 0.5, 2, 10000, 300, 1);
}

void test8()
{
    shock_test(1, 2, 0.1, 1, -2, 0.1,
        25, 0.5, 0.8, 200, 300);
}

void test8a()
{
    shock_test_adaptative(1, 2, 0.1, 1, -2, 0.1,
        25, 0.5, 0.8, 10000, 300, 1);
}

void test9()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 586.6430139510742;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 195243.1507335173;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 1000;
    auto A = [](double x){return 2 + 2./1000*x;};
    auto rho_ini = [](double x){return 2.14;};
    auto v_ini = [=](double x){return 585.1/rho_ini(x)/A(x);};
    auto P_ini = [](double x){return 180000;};
    unsigned int i = 28;
    double t_end = 100;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 0.4, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);

    Simulation sim1("test" + std::to_string(i) + "-HLL", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLL);
    sim1.simulate_until(t_end);
    sim1.write_to_file("test" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("test" + std::to_string(i) + "-HLLC", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLLC);
    sim2.simulate_until(t_end);
    sim2.write_to_file("test" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g21(mesh.T(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g41(mesh.A(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g51(mesh.v(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g61(mesh.M(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    g41.legend = "Godunov-HLL";
    g51.legend = "Godunov-HLL";
    g61.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g22(mesh.T(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g42(mesh.A(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g52(mesh.v(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g62(mesh.M(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    g42.legend = "Godunov-HLLC";
    g52.legend = "Godunov-HLLC";
    g62.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g41, &Y_axis4, &X_axis4);
    G4.add(&g42, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("T" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("T" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g51, &Y_axis5, &X_axis5);
    G5.add(&g52, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("T" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("T" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g61, &Y_axis6, &X_axis6);
    G6.add(&g62, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("T" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("T" + std::to_string(i+5) + "small.pdf");
}

void test9a()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 586.6430139510742;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 195243.1507335173;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 1000;
    auto A = [](double x){return 2 + 2./1000*x;};
    auto rho_ini = [](double x){return 2.14;};
    auto v_ini = [=](double x){return 585.1/rho_ini(x)/A(x);};
    auto P_ini = [](double x){return 180000;};
    unsigned int i = 28;
    double t_end = 100;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 0.4, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);

    Simulation sim1("testa" + std::to_string(i) + "-HLL", QR,
                0, L, N_cells, true, 2, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLL);
    sim1.adaptive_refinement_period = 5000;
    sim1.simulate_until(t_end);
    sim1.write_to_file("testa" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("testa" + std::to_string(i) + "-HLLC", QR,
                0, L, N_cells, true, 2, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLLC);
    sim2.adaptive_refinement_period = 5000;
    sim2.simulate_until(t_end);
    sim2.write_to_file("testa" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, true, 2, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.adaptive_refinement_period = 5000;
    sim3.simulate_until(t_end);
    sim3.write_to_file("testa" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g21(mesh.T(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g41(mesh.A(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g51(mesh.v(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g61(mesh.M(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    g41.legend = "Godunov-HLL";
    g51.legend = "Godunov-HLL";
    g61.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g22(mesh.T(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g42(mesh.A(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g52(mesh.v(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g62(mesh.M(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    g42.legend = "Godunov-HLLC";
    g52.legend = "Godunov-HLLC";
    g62.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g41, &Y_axis4, &X_axis4);
    G4.add(&g42, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("A" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("A" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g51, &Y_axis5, &X_axis5);
    G5.add(&g52, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("A" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("A" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g61, &Y_axis6, &X_axis6);
    G6.add(&g62, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("A" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("A" + std::to_string(i+5) + "small.pdf");

    // Graphic* GMP1 = sim1.mesh_plot();
    // GMP1->render_to_scene().render("M" + std::to_string(i) + "P1.pdf");
    // Graphic* GMP2 = sim2.mesh_plot();
    // GMP2->render_to_scene().render("M" + std::to_string(i) + "P2.pdf");
    // Graphic* GMP3 = sim3.mesh_plot();
    // GMP3->render_to_scene().render("M" + std::to_string(i) + "P3.pdf");
}

void test10()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_DENSITY_SPEED_PRESSURE, GasBoundaryConditionsType::FREE,
            [](double t){return 1.232248582752462;}, [](double t){return 367.3316399837257;}, [](double t){return 82475.40846222841;},
            [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 1000;
    auto A = [](double x){return 2 + 2./1000*x;};
    auto rho_ini = [](double x){return 1.2;};
    auto v_ini = [=](double x){return 500;};
    auto P_ini = [](double x){return 100000;};
    unsigned int i = 34;
    double t_end = 100;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 1.2, 2, QR.gamma, Solvers::Gas::SolutionType::SUPERSONIC);

    Simulation sim1("test" + std::to_string(i) + "-HLL", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLL);
    sim1.simulate_until(t_end);
    sim1.write_to_file("test" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("test" + std::to_string(i) + "-HLLC", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLLC);
    sim2.simulate_until(t_end);
    sim2.write_to_file("test" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;
    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g21(mesh.T(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g41(mesh.A(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g51(mesh.v(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g61(mesh.M(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    g41.legend = "Godunov-HLL";
    g51.legend = "Godunov-HLL";
    g61.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g22(mesh.T(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g42(mesh.A(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g52(mesh.v(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g62(mesh.M(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    g42.legend = "Godunov-HLLC";
    g52.legend = "Godunov-HLLC";
    g62.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g41, &Y_axis4, &X_axis4);
    G4.add(&g42, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("T" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("T" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g51, &Y_axis5, &X_axis5);
    G5.add(&g52, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("T" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("T" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g61, &Y_axis6, &X_axis6);
    G6.add(&g62, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("T" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("T" + std::to_string(i+5) + "small.pdf");
}

void test10a()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_DENSITY_SPEED_PRESSURE, GasBoundaryConditionsType::FREE,
            [](double t){return 1.232248582752462;}, [](double t){return 367.3316399837257;}, [](double t){return 82475.40846222841;},
            [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 1000;
    auto A = [](double x){return 2 + 2./1000*x;};
    auto rho_ini = [](double x){return 1.2;};
    auto v_ini = [=](double x){return 500;};
    auto P_ini = [](double x){return 100000;};
    unsigned int i = 34;
    double t_end = 50;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 1.2, 2, QR.gamma, Solvers::Gas::SolutionType::SUPERSONIC);

    Simulation sim1("testa" + std::to_string(i) + "-HLL", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLL);
    sim1.adaptive_refinement_period = 100;
    sim1.simulate_until(t_end);
    sim1.write_to_file("testa" + std::to_string(i) + "-HLL.sim");
    
    Simulation sim2("testa" + std::to_string(i) + "-HLLC", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::HLLC);
    sim2.adaptive_refinement_period = 100;
    sim2.simulate_until(t_end);
    sim1.write_to_file("testa" + std::to_string(i) + "-HLLC.sim");
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.adaptive_refinement_period = 100;
    sim3.simulate_until(t_end);
    sim1.write_to_file("testa" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;
    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim1.instant_gas_mesh;
    AveragePlot g11(mesh.rho(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g21(mesh.T(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g31(mesh.P(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g41(mesh.A(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g51(mesh.v(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    AveragePlot g61(mesh.M(), mesh.x_partition(), Color::BLUE, 1.2*LineWidth::THIN);
    g11.legend = "Godunov-HLL";
    g21.legend = "Godunov-HLL";
    g31.legend = "Godunov-HLL";
    g41.legend = "Godunov-HLL";
    g51.legend = "Godunov-HLL";
    g61.legend = "Godunov-HLL";
    mesh = sim2.instant_gas_mesh;
    AveragePlot g12(mesh.rho(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g22(mesh.T(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g32(mesh.P(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g42(mesh.A(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g52(mesh.v(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    AveragePlot g62(mesh.M(), mesh.x_partition(), Color::RED, 1.1*LineWidth::THIN);
    g12.legend = "Godunov-HLLC";
    g22.legend = "Godunov-HLLC";
    g32.legend = "Godunov-HLLC";
    g42.legend = "Godunov-HLLC";
    g52.legend = "Godunov-HLLC";
    g62.legend = "Godunov-HLLC";
    mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g11, &Y_axis1, &X_axis1);
    G1.add(&g12, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g21, &Y_axis2, &X_axis2);
    G2.add(&g22, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g31, &Y_axis3, &X_axis3);
    G3.add(&g32, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g41, &Y_axis4, &X_axis4);
    G4.add(&g42, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("A" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("A" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g51, &Y_axis5, &X_axis5);
    G5.add(&g52, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("A" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("A" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g61, &Y_axis6, &X_axis6);
    G6.add(&g62, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("A" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("A" + std::to_string(i+5) + "small.pdf");
}

void test11()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FREE,
            [](double t){return 932.844522204161;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 500;};
    auto P_ini = [](double x){return 100000;};
    unsigned int i = 40;
    double t_end = 0.1;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS_super(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUPERSONIC);
    Solvers::Gas::ExactSteadySolver SS_sub(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        if (x <= 1)
        {
            rho_0.push_back(SS_sub.rho(A(x)));
            v_0.push_back(SS_sub.v(A(x)));
            T_0.push_back(SS_sub.T(A(x)));
            P_0.push_back(SS_sub.P(A(x)));
            M_0.push_back(SS_sub.M(A(x)));
        }
        else
        {
            rho_0.push_back(SS_super.rho(A(x)));
            v_0.push_back(SS_super.v(A(x)));
            T_0.push_back(SS_super.T(A(x)));
            P_0.push_back(SS_super.P(A(x)));
            M_0.push_back(SS_super.M(A(x)));
        }
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("T" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("T" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("T" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("T" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("T" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("T" + std::to_string(i+5) + "small.pdf");
}

void test11a()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FREE,
            [](double t){return 932.844522204161;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 500;};
    auto P_ini = [](double x){return 100000;};
    unsigned int i = 40;
    double t_end = 0.1;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS_super(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUPERSONIC);
    Solvers::Gas::ExactSteadySolver SS_sub(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.adaptive_refinement_period = 100;
    sim3.simulate_until(t_end);
    sim3.write_to_file("testa" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        if (x <= 1)
        {
            rho_0.push_back(SS_sub.rho(A(x)));
            v_0.push_back(SS_sub.v(A(x)));
            T_0.push_back(SS_sub.T(A(x)));
            P_0.push_back(SS_sub.P(A(x)));
            M_0.push_back(SS_sub.M(A(x)));
        }
        else
        {
            rho_0.push_back(SS_super.rho(A(x)));
            v_0.push_back(SS_super.v(A(x)));
            T_0.push_back(SS_super.T(A(x)));
            P_0.push_back(SS_super.P(A(x)));
            M_0.push_back(SS_super.M(A(x)));
        }
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("A" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("A" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("A" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("A" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("A" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("A" + std::to_string(i+5) + "small.pdf");
}

void test12()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 932.845070324472;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 187894.3684832251;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 932.77/2.22/A(x);};
    auto P_ini = [](double x){return 180000;};
    unsigned int i = 46;
    double t_end = 0.2;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;
    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, 2*LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, 2*LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, 2*LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, 2*LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, 2*LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, 2*LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("T" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("T" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("T" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("T" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("T" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("T" + std::to_string(i+5) + "small.pdf");
}

void test12a()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 932.845070324472;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 187894.3684832251;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 100/A(x);};
    auto P_ini = [](double x){return 180000;};
    unsigned int i = 46;
    double t_end = 0.1;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;

    Solvers::Gas::ExactSteadySolver SS(2.32, 200000, 300, 1, 2, QR.gamma, Solvers::Gas::SolutionType::SUBSONIC);
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe,
                [] (std::function<double(double x)> f, const double a, const double b)
                    {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
                [] (const GasCell& C, const double t){return 0.;},
                500, 0.0136, 0.00136, 1./100, 1./2000, 1./1000);
    sim3.adaptive_refinement_period = 1000;
    sim3.simulate_until(t_end);
    sim3.write_to_file("testa" + std::to_string(i) + "-Roe");

    std::vector<double> X;
    std::vector<double> A_0;
    std::vector<double> rho_0;
    std::vector<double> v_0;
    std::vector<double> T_0;
    std::vector<double> P_0;
    std::vector<double> M_0;
    double delta_x = L / (N_points - 1);
    double x = 0;
    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        A_0.push_back(A(x));
        rho_0.push_back(SS.rho(A(x)));
        v_0.push_back(SS.v(A(x)));
        T_0.push_back(SS.T(A(x)));
        P_0.push_back(SS.P(A(x)));
        M_0.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot g10(rho_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g20(T_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g30(P_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g40(A_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g50(v_0, X, Color::BLACK, LineWidth::THIN);
    LinePlot g60(M_0, X, Color::BLACK, LineWidth::THIN);
    g10.legend = "solución exacta";
    g20.legend = "solución exacta";
    g30.legend = "solución exacta";
    g40.legend = "solución exacta";
    g50.legend = "solución exacta";
    g60.legend = "solución exacta";
    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g10, &Y_axis1, &X_axis1);
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g20, &Y_axis2, &X_axis2);
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g30, &Y_axis3, &X_axis3);
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g40, &Y_axis4, &X_axis4);
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("A" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("A" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g50, &Y_axis5, &X_axis5);
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("A" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("A" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g60, &Y_axis6, &X_axis6);
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("A" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("A" + std::to_string(i+5) + "small.pdf");
}

void test13()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 932.845070324472;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 130000;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 932.77/2.22/A(x);};
    auto P_ini = [](double x){return 80000;};
    unsigned int i = 52;
    double t_end = 0.1;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;
    
    Simulation sim3("test" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, false, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.simulate_until(t_end);
    sim3.write_to_file("test" + std::to_string(i) + "-Roe.sim");

    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("T" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("T" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("T" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("T" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("T" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("T" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("T" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("T" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("T" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("T" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("T" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("T" + std::to_string(i+5) + "small.pdf");
}

void test13a()
{
    SolidGasReaction QR(0, 0, 0, 717.5, 287, 0, 101325, 0, 0);
    GasBoundaryConditions BC(GasBoundaryConditionsType::FIXED_FLOW_ENTHALPY, GasBoundaryConditionsType::FIXED_PRESSURE,
            [](double t){return 932.845070324472;}, [](double t){return 301724.1379310346;}, [](double t){return 0;},
            [](double t){return 130000;}, [](double t){return 0;}, [](double t){return 0;});
    double L = 2;
    auto A = [](double x){return 2 + 4*(4.0702 - 2)/4*(x-1)*(x-1);};
    auto rho_ini = [](double x){return 2.22;};
    auto v_ini = [=](double x){return 932.77/2.22/A(x);};
    auto P_ini = [](double x){return 80000;};
    unsigned int i = 52;
    double t_end = 0.1;
    unsigned int N_cells = 1000;
    unsigned int N_points = 300;
    
    Simulation sim3("testa" + std::to_string(i) + "-Roe", QR,
                0, L, N_cells, true, 3, 0.9, 0, BC,
                rho_ini,
                v_ini,
                P_ini,
                A,
                Solvers::Gas::Roe);
    sim3.adaptive_refinement_period = 100;
    sim3.simulate_until(t_end);
    sim3.write_to_file("testa" + std::to_string(i) + "-Roe.sim");

    GasMesh mesh = sim3.instant_gas_mesh;
    AveragePlot g13(mesh.rho(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g23(mesh.T(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g33(mesh.P(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g43(mesh.A(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g53(mesh.v(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g63(mesh.M(), mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g13.legend = "Godunov-Roe";
    g23.legend = "Godunov-Roe";
    g33.legend = "Godunov-Roe";
    g43.legend = "Godunov-Roe";
    g53.legend = "Godunov-Roe";
    g63.legend = "Godunov-Roe";
    
    Axis X_axis1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis1(AxisType::VERTICAL, "$\\rho_g\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g13, &Y_axis1, &X_axis1);
    G1.render_to_scene().render("A" + std::to_string(i) + ".pdf");
    X_axis1.scale = 0.5;
    X_axis1.N_major_ticks = 5;
    Y_axis1.scale = 0.5;
    Y_axis1.N_major_ticks = 5;
    G1.render_to_scene().render("A" + std::to_string(i) + "small.pdf");
    
    Axis X_axis2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis2(AxisType::VERTICAL, "$T$ (K)");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g23, &Y_axis2, &X_axis2);
    G2.render_to_scene().render("A" + std::to_string(i+1) + ".pdf");
    X_axis2.scale = 0.5;
    X_axis2.N_major_ticks = 5;
    Y_axis2.scale = 0.5;
    Y_axis2.N_major_ticks = 5;
    G2.render_to_scene().render("A" + std::to_string(i+1) + "small.pdf");
    
    Axis X_axis3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis3(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g33, &Y_axis3, &X_axis3);
    G3.render_to_scene().render("A" + std::to_string(i+2) + ".pdf");
    X_axis3.scale = 0.5;
    X_axis3.N_major_ticks = 5;
    Y_axis3.scale = 0.5;
    Y_axis3.N_major_ticks = 5;
    G3.render_to_scene().render("A" + std::to_string(i+2) + "small.pdf");

    Axis X_axis4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis4(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g43, &Y_axis4, &X_axis4);
    G4.render_to_scene().render("A" + std::to_string(i+3) + ".pdf");
    X_axis4.scale = 0.5;
    X_axis4.N_major_ticks = 5;
    Y_axis4.scale = 0.5;
    Y_axis4.N_major_ticks = 5;
    G4.render_to_scene().render("A" + std::to_string(i+3) + "small.pdf");

    Axis X_axis5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis5(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g53, &Y_axis5, &X_axis5);
    G5.render_to_scene().render("A" + std::to_string(i+4) + ".pdf");
    X_axis5.scale = 0.5;
    X_axis5.N_major_ticks = 5;
    Y_axis5.scale = 0.5;
    Y_axis5.N_major_ticks = 5;
    G5.render_to_scene().render("A" + std::to_string(i+4) + "small.pdf");

    Axis X_axis6(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g63, &Y_axis6, &X_axis6);
    G6.render_to_scene().render("A" + std::to_string(i+5) + ".pdf");
    X_axis6.scale = 0.5;
    X_axis6.N_major_ticks = 5;
    Y_axis6.scale = 0.5;
    Y_axis6.N_major_ticks = 5;
    G6.render_to_scene().render("A" + std::to_string(i+5) + "small.pdf");
}

void test14()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
            [](double T, double t){return 0.;}, [](double T, double t){return 0.;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [](double x){return (x < 0.5) ? 273+25 : 273+75;};
    unsigned int i = 1;
    double t_end = 10000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testS" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(273+50);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    Y_axis.set_min_value(273+25);
    Y_axis.set_max_value(273+75);
    G.render_to_scene().render("S" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("S" + std::to_string(i) + "small.pdf");
}

void test14a()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
            [](double T, double t){return 0.;}, [](double T, double t){return 0.;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [](double x){return (x < 0.5) ? 273+25 : 273+75;};
    unsigned int i = 1;
    double t_end = 10000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testSA" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testSA" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(273+50);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    Y_axis.set_min_value(273+25);
    Y_axis.set_max_value(273+75);
    G.render_to_scene().render("SA" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("SA" + std::to_string(i) + "small.pdf");
}

void test15()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [](double T, double t){return 273;}, [](double T, double t){return 373;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [](double x){return 300;};
    unsigned int i = 2;
    double t_end = 10000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testS" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(273 + 100/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("S" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("S" + std::to_string(i) + "small.pdf");
}

void test15a()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [](double T, double t){return 273;}, [](double T, double t){return 373;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [](double x){return 300;};
    unsigned int i = 2;
    double t_end = 10000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testSA" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testSA" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(273 + 100/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("SA" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("SA" + std::to_string(i) + "small.pdf");
}

void test16()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [](double T, double t){return 300;}, [](double T, double t){return 300;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300 + 30*sin(M_PI/L*x);};
    unsigned int i = 3;
    double t_end = 600;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testS" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 30*sin(M_PI/L*x)*exp(-M_PI*M_PI/L/L*QR.alpha*t_end));
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("S" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("S" + std::to_string(i) + "small.pdf");
}

void test16a()
{
    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [](double T, double t){return 300;}, [](double T, double t){return 300;});

    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300 + 30*sin(M_PI/L*x);};
    unsigned int i = 3;
    double t_end = 600;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    Simulation sim1("testSA" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_explicit);
    sim1.simulate_until(t_end);

    Simulation sim2("testSA" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 30*sin(M_PI/L*x)*exp(-M_PI*M_PI/L/L*QR.alpha*t_end));
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    SolidMesh mesh = sim1.instant_solid_mesh;
    AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 1.1*LineWidth::THIN);
    g1.legend = "Euler explícito";
    mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("SA" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("SA" + std::to_string(i) + "small.pdf");
}

void test17()
{
    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300;};
    unsigned int i = 4;
    double t_end = 80000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_GRADIENT,
            [](double T, double t){return 300;}, [=](double T, double t){return T/(2*L);});

    // Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
    //     Solvers::Solid::euler_explicit);
    // sim1.simulate_until(t_end);

    Simulation sim2("testS" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 300/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    // SolidMesh mesh = sim1.instant_solid_mesh;
    // AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 2*LineWidth::THIN);
    // g1.legend = "Euler explícito";
    SolidMesh mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    // G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("S" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("S" + std::to_string(i) + "small.pdf");
}

void test17a()
{
    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300;};
    unsigned int i = 4;
    double t_end = 80000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_TEMPERATURE, SolidBoundaryConditionsType::FIXED_GRADIENT,
            [](double T, double t){return 300;}, [=](double T, double t){return T/(2*L);});

    // Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
    //     Solvers::Solid::euler_explicit);
    // sim1.simulate_until(t_end);

    Simulation sim2("testSA" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 300/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    // SolidMesh mesh = sim1.instant_solid_mesh;
    // AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 2*LineWidth::THIN);
    // g1.legend = "Euler explícito";
    SolidMesh mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    // G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("SA" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("SA" + std::to_string(i) + "small.pdf");
}

void test18()
{
    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300;};
    unsigned int i = 5;
    double t_end = 40000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [=](double T, double t){return T/(2*L);}, [](double T, double t){return 450;});

    // Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
    //     Solvers::Solid::euler_explicit);
    // sim1.simulate_until(t_end);

    Simulation sim2("testS" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 150/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    // SolidMesh mesh = sim1.instant_solid_mesh;
    // AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 2*LineWidth::THIN);
    // g1.legend = "Euler explícito";
    SolidMesh mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    // G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("S" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("S" + std::to_string(i) + "small.pdf");
}


void test18a()
{
    double L = 1;
    auto A = [](double x){return 0.01;};
    auto T = [=](double x){return 300;};
    unsigned int i = 5;
    double t_end = 40000;
    unsigned int N_cells = 100;
    unsigned int N_points = 300;

    SolidGasReaction QR(2700, 237, 890, 0, 0,0, 101325, 0, 0);
    SolidBoundaryConditions BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_TEMPERATURE,
            [=](double T, double t){return T/(2*L);}, [](double T, double t){return 450;});

    // Simulation sim1("testS" + std::to_string(i) + "-explicit", QR, 0, L, N_cells, false, 1, 0.9, 0, BC, T, A,
    //     Solvers::Solid::euler_explicit);
    // sim1.simulate_until(t_end);

    Simulation sim2("testSA" + std::to_string(i) + "-implicit", QR, 0, L, N_cells, true, 1, 0.9, 0, BC, T, A,
        Solvers::Solid::euler_implicit);
    for (unsigned int i = 0; (double)i < t_end; i++)
    {
        sim2.update(1);
    }
    
    std::vector<double> X;
    std::vector<double> T_0;
    double delta_x = L / (N_points - 1);
    double x = 0;    
    for (unsigned int i = 0; i < N_points; ++i)
    {
        X.push_back(x);
        T_0.push_back(300 + 150/L*x);
        x += delta_x;
    }
    LinePlot g0(T_0, X, Color::BLACK, LineWidth::THIN);
    g0.legend = "solución exacta";
    // SolidMesh mesh = sim1.instant_solid_mesh;
    // AveragePlot g1(mesh.T(), mesh.x_partition(), Color::BLUE, 2*LineWidth::THIN);
    // g1.legend = "Euler explícito";
    SolidMesh mesh = sim2.instant_solid_mesh;
    AveragePlot g2(mesh.T(), mesh.x_partition(), Color::RED, LineWidth::THIN);
    g2.legend = "Euler implícito";

    Axis X_axis(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis(AxisType::VERTICAL, "$T$ (K)");
    Graphic G;
    G.show_legend = true;
    G.legend_position = LegendPosition::ABOVE;
    G.add(&g0, &Y_axis, &X_axis);
    // G.add(&g1, &Y_axis, &X_axis);
    G.add(&g2, &Y_axis, &X_axis);
    G.render_to_scene().render("SA" + std::to_string(i) + ".pdf");
    X_axis.scale = 0.5;
    X_axis.N_major_ticks = 5;
    Y_axis.scale = 0.5;
    Y_axis.N_major_ticks = 5;
    G.render_to_scene().render("SA" + std::to_string(i) + "small.pdf");
}

void test19()
{
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    double L = 0.3;
    double x_q = 0.29;
    auto v = [](double x){return 0;};
    auto P = [](double x){return 101325;};
    auto T = [](double x){return 300;};
    auto A = [](double x){return M_PI*0.05*0.05;};
    unsigned int i = 1;
    double t_end = 2.15;
    unsigned int N_cells_solid = 300;
    unsigned int N_cells_gas = 10;
    unsigned int N_points = 300;

    SolidBoundaryConditions solid_BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
        [](double T, double t){return 0;}, [](double T, double t){return 0;});
    GasBoundaryConditions gas_BC(GasBoundaryConditionsType::WALL, GasBoundaryConditionsType::WALL);

    Simulation sim3("testB" + std::to_string(i) + "-Roe", QR, 0, L, N_cells_solid, N_cells_gas, true, false, 1, 0.9,
        10000, solid_BC, gas_BC, x_q, v, P, T, A, Solvers::Gas::Roe, Solvers::Solid::euler_implicit,
        [] (std::function<double(double x)> f, const double a, const double b)
            {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        [] (const GasCell& C, const double t){return 0.;},
        10, 0.005, 0.0005, 1./50, 1./10000, 1./10000);
    sim3.simulate_until(t_end);
    sim3.write_to_file("testB" + std::to_string(i) + "-Roe.sim");

    SolidMesh solid_mesh = sim3.instant_solid_mesh;
    GasMesh gas_mesh = sim3.instant_gas_mesh;
    AveragePlot g1sRoe(solid_mesh.T(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g1gRoe(gas_mesh.T(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g2Roe(gas_mesh.v(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g3Roe(gas_mesh.M(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g4Roe(gas_mesh.P(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5Roe(gas_mesh.rho(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g1sRoe.legend = "Roe";
    g2Roe.legend = "Roe";
    g3Roe.legend = "Roe";
    g4Roe.legend = "Roe";
    g5Roe.legend = "Roe";

    Axis X_axis_1(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_1(AxisType::VERTICAL, "$T$ (K)");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g1sRoe, &Y_axis_1, &X_axis_1);
    G1.add(&g1gRoe, &Y_axis_1, &X_axis_1);
    G1.render_to_scene().render("B" + std::to_string(i) + ".pdf");

    Axis X_axis_2(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_2(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g2Roe, &Y_axis_2, &X_axis_2);
    X_axis_2.set_min_value(0);
    G2.render_to_scene().render("B" + std::to_string(i+1) + ".pdf");

    Axis X_axis_3(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_3(AxisType::VERTICAL, "$M$");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g3Roe, &Y_axis_3, &X_axis_3);
    X_axis_3.set_min_value(0);
    G3.render_to_scene().render("B" + std::to_string(i+2) + ".pdf");

    Axis X_axis_4(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_4(AxisType::VERTICAL, "$P$ (Pa)");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g4Roe, &Y_axis_4, &X_axis_4);
    X_axis_4.set_min_value(0);
    G4.render_to_scene().render("B" + std::to_string(i+3) + ".pdf");

    Axis X_axis_5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_5(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g5Roe, &Y_axis_5, &X_axis_5);
    X_axis_5.set_min_value(0);
    G5.render_to_scene().render("B" + std::to_string(i+4) + ".pdf");
}

void test20()
{
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    double L = 0.3;
    double x_q = 0.29;
    auto v = [](double x){return 0;};
    auto P = [](double x){return 101325;};
    auto T = [](double x)
    {
        if (x > 0.285 && x < 0.29)
        {
            return 1000;
        }
        else
        {
            return 300;
        }
    };
    auto A = [](double x){return M_PI*0.05*0.05;};
    unsigned int i = 6;
    double t_end = 35;
    unsigned int N_cells_solid = 10000;
    unsigned int N_cells_gas = 10;
    unsigned int N_points = 300;

    SolidBoundaryConditions solid_BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
        [](double T, double t){return 0;}, [](double T, double t){return 0;});
    GasBoundaryConditions gas_BC(GasBoundaryConditionsType::WALL, GasBoundaryConditionsType::FIXED_PRESSURE,
        [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;},
        [](double t){return 101325;}, [](double t){return 0;}, [](double t){return 0;});

    Simulation sim3("testB" + std::to_string(i) + "-Roe", QR, 0, L, N_cells_solid, N_cells_gas, true, false, 1, 0.9,
        10000, solid_BC, gas_BC, x_q, v, P, T, A, Solvers::Gas::Roe, Solvers::Solid::euler_implicit,
        [] (std::function<double(double x)> f, const double a, const double b)
            {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        [] (const GasCell& C, const double t){return 0.;},
        10, 0.005, 0.0005, 1./50, 1./10000, 1./10000);
    sim3.simulate_until(t_end);
    sim3.write_to_file("testB" + std::to_string(i) + "-Roe.sim");

    SolidMesh solid_mesh = sim3.instant_solid_mesh;
    GasMesh gas_mesh = sim3.instant_gas_mesh;
    AveragePlot g1s(solid_mesh.T(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g1g(gas_mesh.T(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g2(gas_mesh.v(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g3(gas_mesh.P(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g4(gas_mesh.rho(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5(gas_mesh.M(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g1s.legend = "Roe";
    g2.legend = "Roe";
    g3.legend = "Roe";
    g4.legend = "Roe";
    g5.legend = "Roe";

    Axis X_axis_1(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_1(AxisType::VERTICAL, "$T$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g1s, &Y_axis_1, &X_axis_1);
    G1.add(&g1g, &Y_axis_1, &X_axis_1);
    G1.render_to_scene().render("B" + std::to_string(i) + ".pdf");

    Axis X_axis_2(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_2(AxisType::VERTICAL, "$v$");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g2, &Y_axis_2, &X_axis_2);
    X_axis_2.set_min_value(0);
    G2.render_to_scene().render("B" + std::to_string(i+1) + ".pdf");

    Axis X_axis_3(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_3(AxisType::VERTICAL, "$P$");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g3, &Y_axis_3, &X_axis_3);
    X_axis_3.set_min_value(0);
    G3.render_to_scene().render("B" + std::to_string(i+2) + ".pdf");

    Axis X_axis_4(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_4(AxisType::VERTICAL, "$\\rho$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g4, &Y_axis_4, &X_axis_4);
    X_axis_4.set_min_value(0);
    G4.render_to_scene().render("B" + std::to_string(i+3) + ".pdf");

    Axis X_axis_5(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_5(AxisType::VERTICAL, "$M$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g5, &Y_axis_5, &X_axis_5);
    X_axis_5.set_min_value(0);
    G5.render_to_scene().render("B" + std::to_string(i+4) + ".pdf");
}

void test21()
{
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    double L = 0.35;
    double A_c = M_PI*0.05*0.05;
    double A_e = A_c / 10;
    double x_c = 0.3;
    double x_e = 0.33;
    double L_t = x_e - x_c;
    double x_m = (x_c + x_e) / 2;
    double x_q = 0.29;
    auto v = [](double x){return 0;};
    auto P = [](double x){return 101325;};
    auto T = [](double x)
    {
        if (x > 0.285 && x < 0.29)
        {
            return 1000;
        }
        else
        {
            return 300;
        }
    };
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
    unsigned int i = 11;
    double t_end = 30;
    unsigned int N_cells_solid = 300;
    unsigned int N_cells_gas = 50;
    unsigned int N_points = 300;

    SolidBoundaryConditions solid_BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
        [](double T, double t){return 0;}, [](double T, double t){return 0;});
    GasBoundaryConditions gas_BC(GasBoundaryConditionsType::WALL, GasBoundaryConditionsType::FIXED_PRESSURE,
        [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;},
        [](double t){return 101325;}, [](double t){return 0;}, [](double t){return 0;});

    Simulation sim3("testB" + std::to_string(i) + "-Roe", QR, 0, L, N_cells_solid, N_cells_gas, true, false, 1, 0.9,
        10000, solid_BC, gas_BC, x_q, v, P, T, A, Solvers::Gas::Roe, Solvers::Solid::euler_implicit,
        [] (std::function<double(double x)> f, const double a, const double b)
            {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        [] (const GasCell& C, const double t){return 0.;},
        10, 0.005, 0.0005, 1./50, 1./10000, 1./10000);
    sim3.simulate_until(t_end);
    sim3.write_to_file("testB" + std::to_string(i) + "-Roe.sim");

    SolidMesh solid_mesh = sim3.instant_solid_mesh;
    GasMesh gas_mesh = sim3.instant_gas_mesh;
    AveragePlot g1s(solid_mesh.T(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g1g(gas_mesh.T(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g2(gas_mesh.v(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g3(gas_mesh.P(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g4(gas_mesh.rho(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5s(solid_mesh.A(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5g(gas_mesh.A(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g6(gas_mesh.M(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g1s.legend = "Roe";
    g2.legend = "Roe";
    g3.legend = "Roe";
    g4.legend = "Roe";
    g5s.legend = "Roe";
    g6.legend = "Roe";

    Axis X_axis_1(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_1(AxisType::VERTICAL, "$T$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g1s, &Y_axis_1, &X_axis_1);
    G1.add(&g1g, &Y_axis_1, &X_axis_1);
    G1.render_to_scene().render("B" + std::to_string(i) + ".pdf");

    Axis X_axis_2(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_2(AxisType::VERTICAL, "$v$");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g2, &Y_axis_2, &X_axis_2);
    X_axis_2.set_min_value(0);
    G2.render_to_scene().render("B" + std::to_string(i+1) + ".pdf");

    Axis X_axis_3(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_3(AxisType::VERTICAL, "$P$");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g3, &Y_axis_3, &X_axis_3);
    X_axis_3.set_min_value(0);
    G3.render_to_scene().render("B" + std::to_string(i+2) + ".pdf");

    Axis X_axis_4(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_4(AxisType::VERTICAL, "$\\rho$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g4, &Y_axis_4, &X_axis_4);
    X_axis_4.set_min_value(0);
    G4.render_to_scene().render("B" + std::to_string(i+3) + ".pdf");

    Axis X_axis_5(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_5(AxisType::VERTICAL, "$A$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g5s, &Y_axis_5, &X_axis_5);
    G5.add(&g5g, &Y_axis_5, &X_axis_5);
    G5.render_to_scene().render("B" + std::to_string(i+4) + ".pdf");
    
    Axis X_axis_6(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g6, &Y_axis_6, &X_axis_6);
    G6.render_to_scene().render("B" + std::to_string(i+5) + ".pdf");
}

void test22()
{
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    double L = 0.35;
    double A_c = M_PI*0.05*0.05;
    double A_e = A_c * 10;
    double x_c = 0.3;
    double x_e = 0.33;
    double L_t = x_e - x_c;
    double x_m = (x_c + x_e) / 2;
    double x_q = 0.29;
    auto v = [](double x){return 0;};
    auto P = [](double x){return 101325;};
    auto T = [](double x)
    {
        if (x > 0.285 && x < 0.29)
        {
            return 1000;
        }
        else
        {
            return 300;
        }
    };
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
    unsigned int i = 17;
    double t_end = 35;
    unsigned int N_cells_solid = 300;
    unsigned int N_cells_gas = 50;
    unsigned int N_points = 300;

    SolidBoundaryConditions solid_BC(SolidBoundaryConditionsType::FIXED_GRADIENT, SolidBoundaryConditionsType::FIXED_GRADIENT,
        [](double T, double t){return 0;}, [](double T, double t){return 0;});
    GasBoundaryConditions gas_BC(GasBoundaryConditionsType::WALL, GasBoundaryConditionsType::FIXED_PRESSURE,
        [](double t){return 0;}, [](double t){return 0;}, [](double t){return 0;},
        [](double t){return 101325;}, [](double t){return 0;}, [](double t){return 0;});

    Simulation sim3("testB" + std::to_string(i) + "-Roe", QR, 0, L, N_cells_solid, N_cells_gas, true, false, 1, 0.9,
        10000, solid_BC, gas_BC, x_q, v, P, T, A, Solvers::Gas::Roe, Solvers::Solid::euler_implicit,
        [] (std::function<double(double x)> f, const double a, const double b)
            {return Math::Integrators::Gauss_Konrad_G7_K15(f, a, b);},
        [] (const GasCell& C, const double t){return 0.;},
        10, 0.005, 0.0005, 1./50, 1./10000, 1./10000);
    sim3.simulate_until(t_end);
    sim3.write_to_file("testB" + std::to_string(i) + "-Roe.sim");

    SolidMesh solid_mesh = sim3.instant_solid_mesh;
    GasMesh gas_mesh = sim3.instant_gas_mesh;
    AveragePlot g1s(solid_mesh.T(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g1g(gas_mesh.T(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g2(gas_mesh.v(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g3(gas_mesh.P(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g4(gas_mesh.rho(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5s(solid_mesh.A(), solid_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g5g(gas_mesh.A(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    AveragePlot g6(gas_mesh.M(), gas_mesh.x_partition(), Color::GREEN, LineWidth::THIN);
    g1s.legend = "Roe";
    g2.legend = "Roe";
    g3.legend = "Roe";
    g4.legend = "Roe";
    g5s.legend = "Roe";
    g6.legend = "Roe";

    Axis X_axis_1(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_1(AxisType::VERTICAL, "$T$");
    Graphic G1;
    G1.show_legend = true;
    G1.legend_position = LegendPosition::ABOVE;
    G1.add(&g1s, &Y_axis_1, &X_axis_1);
    G1.add(&g1g, &Y_axis_1, &X_axis_1);
    G1.render_to_scene().render("B" + std::to_string(i) + ".pdf");

    Axis X_axis_2(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_2(AxisType::VERTICAL, "$v$");
    Graphic G2;
    G2.show_legend = true;
    G2.legend_position = LegendPosition::ABOVE;
    G2.add(&g2, &Y_axis_2, &X_axis_2);
    X_axis_2.set_min_value(0);
    G2.render_to_scene().render("B" + std::to_string(i+1) + ".pdf");

    Axis X_axis_3(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_3(AxisType::VERTICAL, "$P$");
    Graphic G3;
    G3.show_legend = true;
    G3.legend_position = LegendPosition::ABOVE;
    G3.add(&g3, &Y_axis_3, &X_axis_3);
    X_axis_3.set_min_value(0);
    G3.render_to_scene().render("B" + std::to_string(i+2) + ".pdf");

    Axis X_axis_4(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_4(AxisType::VERTICAL, "$\\rho$");
    Graphic G4;
    G4.show_legend = true;
    G4.legend_position = LegendPosition::ABOVE;
    G4.add(&g4, &Y_axis_4, &X_axis_4);
    X_axis_4.set_min_value(0);
    G4.render_to_scene().render("B" + std::to_string(i+3) + ".pdf");

    Axis X_axis_5(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_5(AxisType::VERTICAL, "$A$");
    Graphic G5;
    G5.show_legend = true;
    G5.legend_position = LegendPosition::ABOVE;
    G5.add(&g5s, &Y_axis_5, &X_axis_5);
    G5.add(&g5g, &Y_axis_5, &X_axis_5);
    G5.render_to_scene().render("B" + std::to_string(i+4) + ".pdf");

    Axis X_axis_6(AxisType::HORIZONTAL, "$x$");
    Axis Y_axis_6(AxisType::VERTICAL, "$M$");
    Graphic G6;
    G6.show_legend = true;
    G6.legend_position = LegendPosition::ABOVE;
    G6.add(&g6, &Y_axis_6, &X_axis_6);
    G6.render_to_scene().render("B" + std::to_string(i+5) + ".pdf");
}