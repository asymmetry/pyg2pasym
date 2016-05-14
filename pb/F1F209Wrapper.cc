// -*- C++ -*-

/* class F1F209Wrapper
 * Wrapper class to get f1, f2, and xs with P. Bosted fitting.
 * Unit of xs is ub/MeV-sr.
 * Valid for all W<5 GeV and all Q2<11 GeV2.
 * Use -fno-leading-underscore -fno-second-underscore when compiling F1F209.f
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   May 2014, Comments modified By Jixie Zhang

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "F1F209Wrapper.hh"

using namespace std;

extern "C" {
void f1f2in09_(double* Z, double* A, double* Q2, double* W2, double* F1, double* F2, double* rc);
void f1f2qe09_(double* Z, double* A, double* Q2, double* W2, double* F1, double* F2);
}

F1F209Wrapper::F1F209Wrapper()
{
    // Nothing to do
}

F1F209Wrapper::~F1F209Wrapper()
{
    // Nothing to do
}

void F1F209Wrapper::GetF1F2IN09(double Z, double A, double Q2, double W2, double &F1, double &F2, double &rc)
{
    double F1o, F2o, rco;
    f1f2in09_(&Z, &A, &Q2, &W2, &F1o, &F2o, &rco);
    F1 = F1o;
    F2 = F2o;
    rc = rco;
}

void F1F209Wrapper::GetF1F2QE09(double Z, double A, double Q2, double W2, double &F1, double &F2)
{
    double F1o, F2o, rco;
    f1f2qe09_(&Z, &A, &Q2, &W2, &F1o, &F2o);
    F1 = F1o;
    F2 = F2o;
}

double F1F209Wrapper::GetXS(double Z, double A, double Ei, double Ef, double theta)
{
    const double M = 0.93825;

    double nu = Ei - Ef;
    double Q2 = 4. * Ei * Ef * (sin(abs(theta) / 2.))*(sin(abs(theta) / 2.));
    double w2 = M * M + 2. * M * nu - Q2;

    double F1, F2, r;
    double xs1, xs2;

    GetF1F2IN09(Z, A, Q2, w2, F1, F2, r);

    xs1 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott
    xs1 = xs1 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
    xs1 = xs1 * 389.379;

    GetF1F2QE09(Z, A, Q2, w2, F1, F2);
    xs2 = (2. / 137. * Ef / Q2 * cos(theta / 2.))*(2. / 137. * Ef / Q2 * cos(theta / 2.)); // mott
    xs2 = xs2 * (2. / M * F1 * tan(abs(theta) / 2.) * tan(abs(theta) / 2.) + F2 / nu);
    xs2 = xs2 * 389.379;

    return (xs1 + xs2) / 1000.; // ub/MeV-sr
}
