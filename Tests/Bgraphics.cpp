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
using namespace CPGF::AffineSpace;
using namespace CPGF;

void graphic19()
{
    Simulation sim("testB1-Roe.sim");

    Graphic T;
    Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
    T.show_legend = true;
    T.legend_position = LegendPosition::ABOVE;
    Graphic P;
    Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
    P.show_legend = true;
    P.legend_position = LegendPosition::ABOVE;
    Graphic rho;
    Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    rho.show_legend = true;
    rho.legend_position = LegendPosition::ABOVE;
    Graphic v;
    Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Y_axis_v.set_min_value(-1);
    Y_axis_v.set_max_value(1);
    v.show_legend = true;
    v.legend_position = LegendPosition::ABOVE;

    LinePlot gT_approx(std::vector<Point2d>({Point2d(0, 1735.29), Point2d(0.3, 1735.29)}), Color::BLACK);
    T.add(&gT_approx, &Y_axis_T, &X_axis_T);
    gT_approx.legend = "aproximación estacionaria";
    LinePlot gP_approx(std::vector<Point2d>({Point2d(0, 558.79e6), Point2d(0.3, 558.79e6)}), Color::BLACK);
    P.add(&gP_approx, &Y_axis_P, &X_axis_P);
    gP_approx.legend = "aproximación estacionaria";
    LinePlot grho_approx(std::vector<Point2d>({Point2d(0, 2012.6), Point2d(0.3, 2012.6)}), Color::BLACK);
    rho.add(&grho_approx, &Y_axis_rho, &X_axis_rho);
    grho_approx.legend = "aproximación estacionaria";
    LinePlot gv_approx(std::vector<Point2d>({Point2d(0, 0), Point2d(0.3, 0)}), Color::BLACK);
    v.add(&gv_approx, &Y_axis_v, &X_axis_v);
    gv_approx.legend = "aproximación estacionaria";

    for (unsigned int i = 0; i < sim.gas_mesh->size(); i+= 1900)
    {
        GasMesh gas_mesh = (*sim.gas_mesh)[i];
        SolidMesh solid_mesh = (*sim.solid_mesh)[i];

        // // We fix QR.
        // for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
        // {
        //     C->QR = &sim.QR;
        //     C->update();
        // }

        AveragePlot* gT_s = new AveragePlot(solid_mesh.T(), solid_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        AveragePlot* gT_g = new AveragePlot(gas_mesh.T(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        gT_g->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        T.add(gT_s, &Y_axis_T, &X_axis_T);
        T.add(gT_g, &Y_axis_T, &X_axis_T);
        AveragePlot* gP = new AveragePlot(gas_mesh.P(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        P.add(gP, &Y_axis_P, &X_axis_P);
        gP->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* grho = new AveragePlot(gas_mesh.rho(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        rho.add(grho, &Y_axis_rho, &X_axis_rho);
        grho->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gv = new AveragePlot(gas_mesh.v(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        v.add(gv, &Y_axis_v, &X_axis_v);
        gv->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
    }

    T.render_to_scene().render("B1.pdf");
    X_axis_T.scale = 0.5;
    X_axis_T.N_major_ticks = 5;
    Y_axis_T.scale = 0.5;
    Y_axis_T.N_major_ticks = 5;
    T.render_to_scene().render("B1small.pdf");
    P.render_to_scene().render("B2.pdf");
    X_axis_P.scale = 0.5;
    X_axis_P.N_major_ticks = 5;
    Y_axis_P.scale = 0.5;
    Y_axis_P.N_major_ticks = 5;
    P.render_to_scene().render("B2small.pdf");
    rho.render_to_scene().render("B3.pdf");
    X_axis_rho.scale = 0.5;
    X_axis_rho.N_major_ticks = 5;
    Y_axis_rho.scale = 0.5;
    Y_axis_rho.N_major_ticks = 5;
    rho.render_to_scene().render("B3small.pdf");
    v.render_to_scene().render("B4.pdf");
    X_axis_v.scale = 0.5;
    X_axis_v.N_major_ticks = 5;
    Y_axis_v.scale = 0.5;
    Y_axis_v.N_major_ticks = 5;
    v.render_to_scene().render("B4small.pdf");

    Graphic vq;
    Axis X_axis_vq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_vq(AxisType::VERTICAL, "$v_q\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    LinePlot gvq(sim.v_q_array, sim.t_array);
    vq.add(&gvq, &Y_axis_vq, &X_axis_vq);
    vq.render_to_scene().render("B5.pdf");
    X_axis_vq.scale = 0.5;
    X_axis_vq.N_major_ticks = 5;
    Y_axis_vq.scale = 0.5;
    Y_axis_vq.N_major_ticks = 5;
    vq.render_to_scene().render("B5small.pdf");
    Graphic xq;
    Axis X_axis_xq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_xq(AxisType::VERTICAL, "$x_q$ (m)");
    LinePlot gxq(sim.x_q_array, sim.t_array);
    xq.add(&gxq, &Y_axis_xq, &X_axis_xq);
    xq.render_to_scene().render("B6.pdf");
    X_axis_xq.scale = 0.5;
    X_axis_xq.N_major_ticks = 5;
    Y_axis_xq.scale = 0.5;
    Y_axis_xq.N_major_ticks = 5;
    xq.render_to_scene().render("B6small.pdf");
}

void graphic20()
{
    Simulation sim("testB6-Roe.sim");
    SolidGasReaction QR(2082, 82.02, 895, 731, 160, 0.00338, 6894.76, 0.325, 1e6);

    Graphic T;
    Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
    T.show_legend = true;
    T.legend_position = LegendPosition::ABOVE;
    Graphic P;
    Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
    Y_axis_P.set_min_value(0);
    P.show_legend = true;
    P.legend_position = LegendPosition::ABOVE;
    Graphic rho;
    Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    Y_axis_rho.set_min_value(0);
    rho.show_legend = true;
    rho.legend_position = LegendPosition::ABOVE;
    Graphic v;
    Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Y_axis_v.set_min_value(0);
    v.show_legend = true;
    v.legend_position = LegendPosition::ABOVE;

    LinePlot gT_approx(std::vector<Point2d>({Point2d(0, 1423.68), Point2d(0.3, 1423.68)}), Color::BLACK);
    T.add(&gT_approx, &Y_axis_T, &X_axis_T);
    gT_approx.legend = "aproximación estacionaria";
    LinePlot gP_approx(std::vector<Point2d>({Point2d(0, 101325), Point2d(0.3, 101325)}), Color::BLACK);
    P.add(&gP_approx, &Y_axis_P, &X_axis_P);
    gP_approx.legend = "aproximación estacionaria";
    LinePlot grho_approx(std::vector<Point2d>({Point2d(0, 0.44820), Point2d(0.3, 0.44820)}), Color::BLACK);
    rho.add(&grho_approx, &Y_axis_rho, &X_axis_rho);
    grho_approx.legend = "aproximación estacionaria";
    LinePlot gv_approx(std::vector<Point2d>({Point2d(0, 37.892452), Point2d(0.3, 37.892452)}), Color::BLACK);
    v.add(&gv_approx, &Y_axis_v, &X_axis_v);
    gv_approx.legend = "aproximación estacionaria";

    for (unsigned int i = 100; i < sim.gas_mesh->size(); i+= 930)
    {
        GasMesh gas_mesh = (*sim.gas_mesh)[i];
        SolidMesh solid_mesh = (*sim.solid_mesh)[i];

        // We fix QR.
        for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
        {
            C->QR = &QR;
            C->update();
        }

        AveragePlot* gT_s = new AveragePlot(solid_mesh.T(), solid_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        AveragePlot* gT_g = new AveragePlot(gas_mesh.T(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        gT_g->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        T.add(gT_s, &Y_axis_T, &X_axis_T);
        T.add(gT_g, &Y_axis_T, &X_axis_T);
        AveragePlot* gP = new AveragePlot(gas_mesh.P(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        P.add(gP, &Y_axis_P, &X_axis_P);
        gP->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* grho = new AveragePlot(gas_mesh.rho(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        rho.add(grho, &Y_axis_rho, &X_axis_rho);
        grho->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gv = new AveragePlot(gas_mesh.v(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        v.add(gv, &Y_axis_v, &X_axis_v);
        gv->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
    }

    T.render_to_scene().render("B7.pdf");
    X_axis_T.scale = 0.5;
    X_axis_T.N_major_ticks = 5;
    Y_axis_T.scale = 0.5;
    Y_axis_T.N_major_ticks = 5;
    T.render_to_scene().render("B7small.pdf");
    P.render_to_scene().render("B8.pdf");
    X_axis_P.scale = 0.5;
    X_axis_P.N_major_ticks = 5;
    Y_axis_P.scale = 0.5;
    Y_axis_P.N_major_ticks = 5;
    P.render_to_scene().render("B8small.pdf");
    rho.render_to_scene().render("B9.pdf");
    X_axis_rho.scale = 0.5;
    X_axis_rho.N_major_ticks = 5;
    Y_axis_rho.scale = 0.5;
    Y_axis_rho.N_major_ticks = 5;
    rho.render_to_scene().render("B9small.pdf");
    v.render_to_scene().render("B10.pdf");
    X_axis_v.scale = 0.5;
    X_axis_v.N_major_ticks = 5;
    Y_axis_v.scale = 0.5;
    Y_axis_v.N_major_ticks = 5;
    v.render_to_scene().render("B10small.pdf");

    Graphic vq;
    Axis X_axis_vq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_vq(AxisType::VERTICAL, "$v_q\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Y_axis_vq.set_min_value(0);
    LinePlot gvq(sim.v_q_array, sim.t_array);
    vq.add(&gvq, &Y_axis_vq, &X_axis_vq);
    vq.render_to_scene().render("B11.pdf");
    X_axis_vq.scale = 0.5;
    X_axis_vq.N_major_ticks = 5;
    Y_axis_vq.scale = 0.5;
    Y_axis_vq.N_major_ticks = 5;
    vq.render_to_scene().render("B11small.pdf");
    Graphic xq;
    Axis X_axis_xq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_xq(AxisType::VERTICAL, "$x_q$ (m)");
    LinePlot gxq(sim.x_q_array, sim.t_array);
    xq.add(&gxq, &Y_axis_xq, &X_axis_xq);
    xq.render_to_scene().render("B12.pdf");
    X_axis_xq.scale = 0.5;
    X_axis_xq.N_major_ticks = 5;
    Y_axis_xq.scale = 0.5;
    Y_axis_xq.N_major_ticks = 5;
    xq.render_to_scene().render("B12small.pdf");
}

void graphic21()
{
    Simulation sim("testB11-Roe.sim");

    Graphic T;
    Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
    T.show_legend = true;
    T.legend_position = LegendPosition::ABOVE;
    Graphic P;
    Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
    P.show_legend = true;
    P.legend_position = LegendPosition::ABOVE;
    Graphic rho;
    Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    rho.show_legend = true;
    rho.legend_position = LegendPosition::ABOVE;
    Graphic v;
    Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    v.show_legend = true;
    v.legend_position = LegendPosition::ABOVE;
    Graphic M;
    Axis X_axis_M(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_M(AxisType::VERTICAL, "$M$");
    M.show_legend = true;
    M.legend_position = LegendPosition::ABOVE;
    Graphic AG;
    Axis X_axis_AG(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_AG(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Y_axis_AG.set_min_value(0);

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
    std::vector<double> X_approx;
    std::vector<double> A_approx;
    std::vector<double> T_approx;
    std::vector<double> P_approx;
    std::vector<double> rho_approx;
    std::vector<double> v_approx;
    std::vector<double> M_approx;
    unsigned int N_points = 300;
    double x = 0;
    double delta_x = 0.35 / (N_points - 1);
    for (unsigned int i = 0; i < 300; i++)
    {
        X_approx.push_back(x);
        A_approx.push_back(A(x));
        T_approx.push_back(SS.T(A(x)));
        P_approx.push_back(SS.P(A(x)));
        rho_approx.push_back(SS.rho(A(x)));
        v_approx.push_back(SS.v(A(x)));
        M_approx.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot gT_approx(T_approx, X_approx, Color::BLACK);
    T.add(&gT_approx, &Y_axis_T, &X_axis_T);
    gT_approx.legend = "aproximación estacionaria";
    LinePlot gP_approx(P_approx, X_approx, Color::BLACK);
    P.add(&gP_approx, &Y_axis_P, &X_axis_P);
    gP_approx.legend = "aproximación estacionaria";
    LinePlot grho_approx(rho_approx, X_approx, Color::BLACK);
    rho.add(&grho_approx, &Y_axis_rho, &X_axis_rho);
    grho_approx.legend = "aproximación estacionaria";
    LinePlot gv_approx(v_approx, X_approx, Color::BLACK);
    v.add(&gv_approx, &Y_axis_v, &X_axis_v);
    gv_approx.legend = "aproximación estacionaria";
    LinePlot gM_approx(M_approx, X_approx, Color::BLACK);
    M.add(&gM_approx, &Y_axis_M, &X_axis_M);
    gM_approx.legend = "aproximación estacionaria";
    LinePlot gAG_approx(A_approx, X_approx, Color::BLUE);
    AG.add(&gAG_approx, &Y_axis_AG, &X_axis_AG);

    for (unsigned int i = 100; i < sim.gas_mesh->size(); i+= 900)
    {
        GasMesh gas_mesh = (*sim.gas_mesh)[i];
        SolidMesh solid_mesh = (*sim.solid_mesh)[i];

        // // We fix QR.
        // for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
        // {
        //     C->QR = &sim.QR;
        //     C->update();
        // }

        AveragePlot* gT_s = new AveragePlot(solid_mesh.T(), solid_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        AveragePlot* gT_g = new AveragePlot(gas_mesh.T(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        gT_g->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        T.add(gT_s, &Y_axis_T, &X_axis_T);
        T.add(gT_g, &Y_axis_T, &X_axis_T);
        AveragePlot* gP = new AveragePlot(gas_mesh.P(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        P.add(gP, &Y_axis_P, &X_axis_P);
        gP->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* grho = new AveragePlot(gas_mesh.rho(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        rho.add(grho, &Y_axis_rho, &X_axis_rho);
        grho->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gv = new AveragePlot(gas_mesh.v(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        v.add(gv, &Y_axis_v, &X_axis_v);
        gv->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gM = new AveragePlot(gas_mesh.M(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        M.add(gM, &Y_axis_M, &X_axis_M);
        gM->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
    }

    T.render_to_scene().render("B13.pdf");
    X_axis_T.scale = 0.5;
    X_axis_T.N_major_ticks = 5;
    Y_axis_T.scale = 0.5;
    Y_axis_T.N_major_ticks = 5;
    T.render_to_scene().render("B13small.pdf");
    P.render_to_scene().render("B14.pdf");
    X_axis_P.scale = 0.5;
    X_axis_P.N_major_ticks = 5;
    Y_axis_P.scale = 0.5;
    Y_axis_P.N_major_ticks = 5;
    P.render_to_scene().render("B14small.pdf");
    rho.render_to_scene().render("B15.pdf");
    X_axis_rho.scale = 0.5;
    X_axis_rho.N_major_ticks = 5;
    Y_axis_rho.scale = 0.5;
    Y_axis_rho.N_major_ticks = 5;
    rho.render_to_scene().render("B15small.pdf");
    v.render_to_scene().render("B16.pdf");
    X_axis_v.scale = 0.5;
    X_axis_v.N_major_ticks = 5;
    Y_axis_v.scale = 0.5;
    Y_axis_v.N_major_ticks = 5;
    v.render_to_scene().render("B16small.pdf");
    M.render_to_scene().render("B17.pdf");
    X_axis_M.scale = 0.5;
    X_axis_M.N_major_ticks = 5;
    Y_axis_M.scale = 0.5;
    Y_axis_M.N_major_ticks = 5;
    M.render_to_scene().render("B17small.pdf");
    AG.render_to_scene().render("B18.pdf");
    X_axis_AG.scale = 0.5;
    X_axis_AG.N_major_ticks = 5;
    Y_axis_AG.scale = 0.5;
    Y_axis_AG.N_major_ticks = 5;
    AG.render_to_scene().render("B18small.pdf");

    Graphic vq;
    Axis X_axis_vq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_vq(AxisType::VERTICAL, "$v_q\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Y_axis_vq.set_min_value(0);
    LinePlot gvq(sim.v_q_array, sim.t_array);
    vq.add(&gvq, &Y_axis_vq, &X_axis_vq);
    vq.render_to_scene().render("B19.pdf");
    X_axis_vq.scale = 0.5;
    X_axis_vq.N_major_ticks = 5;
    Y_axis_vq.scale = 0.5;
    Y_axis_vq.N_major_ticks = 5;
    vq.render_to_scene().render("B19small.pdf");
    Graphic xq;
    Axis X_axis_xq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_xq(AxisType::VERTICAL, "$x_q$ (m)");
    LinePlot gxq(sim.x_q_array, sim.t_array);
    xq.add(&gxq, &Y_axis_xq, &X_axis_xq);
    xq.render_to_scene().render("B20.pdf");
    X_axis_xq.scale = 0.5;
    X_axis_xq.N_major_ticks = 5;
    Y_axis_xq.scale = 0.5;
    Y_axis_xq.N_major_ticks = 5;
    xq.render_to_scene().render("B20small.pdf");
}

void graphic22()
{
    Simulation sim("testB16-Roe.sim");

    Graphic T;
    Axis X_axis_T(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_T(AxisType::VERTICAL, "$T$ (K)");
    T.show_legend = true;
    T.legend_position = LegendPosition::ABOVE;
    Graphic P;
    Axis X_axis_P(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_P(AxisType::VERTICAL, "$P$ (Pa)");
    P.show_legend = true;
    P.legend_position = LegendPosition::ABOVE;
    Graphic rho;
    Axis X_axis_rho(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_rho(AxisType::VERTICAL, "$\\rho\\:\\left(\\frac{\\text{kg}}{\\text{m}^3}\\right)$");
    rho.show_legend = true;
    rho.legend_position = LegendPosition::ABOVE;
    Graphic v;
    Axis X_axis_v(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_v(AxisType::VERTICAL, "$v\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    v.show_legend = true;
    v.legend_position = LegendPosition::ABOVE;
    Graphic M;
    Axis X_axis_M(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_M(AxisType::VERTICAL, "$M$");
    M.show_legend = true;
    M.legend_position = LegendPosition::ABOVE;
    Graphic AG;
    Axis X_axis_AG(AxisType::HORIZONTAL, "$x$ (m)");
    Axis Y_axis_AG(AxisType::VERTICAL, "$A\\:\\left(\\text{m}^2\\right)$");
    Y_axis_AG.set_min_value(0);

    double L = 0.35;
    double A_c = M_PI*0.05*0.05;
    double A_e = A_c * 10;
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
    Solvers::Gas::ExactSteadySolver SS(0.443025, 101328, 1424.49, 0.072065, A_c, 1.218878, Solvers::Gas::SolutionType::SUBSONIC);
    std::vector<double> X_approx;
    std::vector<double> A_approx;
    std::vector<double> T_approx;
    std::vector<double> P_approx;
    std::vector<double> rho_approx;
    std::vector<double> v_approx;
    std::vector<double> M_approx;
    unsigned int N_points = 300;
    double x = 0;
    double delta_x = 0.35 / (N_points - 1);
    for (unsigned int i = 0; i < 300; i++)
    {
        X_approx.push_back(x);
        A_approx.push_back(A(x));
        T_approx.push_back(SS.T(A(x)));
        P_approx.push_back(SS.P(A(x)));
        rho_approx.push_back(SS.rho(A(x)));
        v_approx.push_back(SS.v(A(x)));
        M_approx.push_back(SS.M(A(x)));
        x += delta_x;
    }
    LinePlot gT_approx(T_approx, X_approx, Color::BLACK);
    T.add(&gT_approx, &Y_axis_T, &X_axis_T);
    gT_approx.legend = "aproximación estacionaria";
    LinePlot gP_approx(P_approx, X_approx, Color::BLACK);
    P.add(&gP_approx, &Y_axis_P, &X_axis_P);
    gP_approx.legend = "aproximación estacionaria";
    LinePlot grho_approx(rho_approx, X_approx, Color::BLACK);
    rho.add(&grho_approx, &Y_axis_rho, &X_axis_rho);
    grho_approx.legend = "aproximación estacionaria";
    LinePlot gv_approx(v_approx, X_approx, Color::BLACK);
    v.add(&gv_approx, &Y_axis_v, &X_axis_v);
    gv_approx.legend = "aproximación estacionaria";
    LinePlot gM_approx(M_approx, X_approx, Color::BLACK);
    M.add(&gM_approx, &Y_axis_M, &X_axis_M);
    gM_approx.legend = "aproximación estacionaria";
    LinePlot gAG_approx(A_approx, X_approx, Color::BLUE);
    AG.add(&gAG_approx, &Y_axis_AG, &X_axis_AG);

    for (unsigned int i = 100; i < sim.gas_mesh->size(); i+= 900)
    {
        GasMesh gas_mesh = (*sim.gas_mesh)[i];
        SolidMesh solid_mesh = (*sim.solid_mesh)[i];

        // // We fix QR.
        // for (GasMesh::Iterator C = gas_mesh.begin(); C != gas_mesh.end(); ++C)
        // {
        //     C->QR = &sim.QR;
        //     C->update();
        // }

        AveragePlot* gT_s = new AveragePlot(solid_mesh.T(), solid_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        AveragePlot* gT_g = new AveragePlot(gas_mesh.T(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        gT_g->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        T.add(gT_s, &Y_axis_T, &X_axis_T);
        T.add(gT_g, &Y_axis_T, &X_axis_T);
        AveragePlot* gP = new AveragePlot(gas_mesh.P(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        P.add(gP, &Y_axis_P, &X_axis_P);
        gP->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* grho = new AveragePlot(gas_mesh.rho(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        rho.add(grho, &Y_axis_rho, &X_axis_rho);
        grho->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gv = new AveragePlot(gas_mesh.v(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        v.add(gv, &Y_axis_v, &X_axis_v);
        gv->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
        AveragePlot* gM = new AveragePlot(gas_mesh.M(), gas_mesh.x_partition(), Color::mix(Color::BLACK, Color::LAWN_GREEN, i/20000.));
        M.add(gM, &Y_axis_M, &X_axis_M);
        gM->legend = "$t = " + Utilities::format_number(sim.t_array[i]) + "$ s";
    }

    T.render_to_scene().render("B21.pdf");
    X_axis_T.scale = 0.5;
    X_axis_T.N_major_ticks = 5;
    Y_axis_T.scale = 0.5;
    Y_axis_T.N_major_ticks = 5;
    T.render_to_scene().render("B21small.pdf");
    P.render_to_scene().render("B22.pdf");
    X_axis_P.scale = 0.5;
    X_axis_P.N_major_ticks = 5;
    Y_axis_P.scale = 0.5;
    Y_axis_P.N_major_ticks = 5;
    P.render_to_scene().render("B22small.pdf");
    rho.render_to_scene().render("B23.pdf");
    X_axis_rho.scale = 0.5;
    X_axis_rho.N_major_ticks = 5;
    Y_axis_rho.scale = 0.5;
    Y_axis_rho.N_major_ticks = 5;
    rho.render_to_scene().render("B23small.pdf");
    v.render_to_scene().render("B24.pdf");
    X_axis_v.scale = 0.5;
    X_axis_v.N_major_ticks = 5;
    Y_axis_v.scale = 0.5;
    Y_axis_v.N_major_ticks = 5;
    v.render_to_scene().render("B24small.pdf");
    M.render_to_scene().render("B25.pdf");
    X_axis_M.scale = 0.5;
    X_axis_M.N_major_ticks = 5;
    Y_axis_M.scale = 0.5;
    Y_axis_M.N_major_ticks = 5;
    M.render_to_scene().render("B25small.pdf");
    AG.render_to_scene().render("B26.pdf");
    X_axis_AG.scale = 0.5;
    X_axis_AG.N_major_ticks = 5;
    Y_axis_AG.scale = 0.5;
    Y_axis_AG.N_major_ticks = 5;
    AG.render_to_scene().render("B26small.pdf");

    Graphic vq;
    Axis X_axis_vq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_vq(AxisType::VERTICAL, "$v_q\\:\\left(\\frac{\\text{m}}{\\text{s}}\\right)$");
    Y_axis_vq.set_min_value(0);
    LinePlot gvq(sim.v_q_array, sim.t_array);
    vq.add(&gvq, &Y_axis_vq, &X_axis_vq);
    vq.render_to_scene().render("B27.pdf");
    X_axis_vq.scale = 0.5;
    X_axis_vq.N_major_ticks = 5;
    Y_axis_vq.scale = 0.5;
    Y_axis_vq.N_major_ticks = 5;
    vq.render_to_scene().render("B27small.pdf");
    Graphic xq;
    Axis X_axis_xq(AxisType::HORIZONTAL, "$t$ (s)");
    Axis Y_axis_xq(AxisType::VERTICAL, "$x_q$ (m)");
    LinePlot gxq(sim.x_q_array, sim.t_array);
    xq.add(&gxq, &Y_axis_xq, &X_axis_xq);
    xq.render_to_scene().render("B28.pdf");
    X_axis_xq.scale = 0.5;
    X_axis_xq.N_major_ticks = 5;
    Y_axis_xq.scale = 0.5;
    Y_axis_xq.N_major_ticks = 5;
    xq.render_to_scene().render("B28small.pdf");
}