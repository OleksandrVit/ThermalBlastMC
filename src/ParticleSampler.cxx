// ROOT
#include <TF1.h>
// Boost
#include <boost/math/quadrature/tanh_sinh.hpp>
// Project
#include "Constants.h"
#include "ParticleSampler.h"

ParticleSampler::ParticleSampler(const smash::ParticleType &p, const double t, const double mu, const double volume, const int N /*= 1000*/, const int Nm /*= 1000*/) : fParticle{p}, fT{t}, fmu{mu}
{
    fUni = std::uniform_real_distribution<double>(0.0, 1.0);
    fa = (p.spin() & 1 ? 1 : -1);
    fmmin = p.min_mass_spectral();
    if (fSF == SpectralFunction::None || p.is_stable())
    {
        fMultiplicity = volume * 0.5 * p.spin_degeneracy() * Density(p.mass()) / (constants::piSqr * constants::hbar3);
        fumin = 0.0;
        fdu = 0.0;
        fdm = 1.0;
        AddCoefficients(fmmin, fmmin + 0.5 * fdm);
    }
    else
    {
        fumin = std::atan(2.0 * (fmmin - p.mass()) / p.width_at_pole());
        const double umax = constants::halfPi;
        const double mmax = p.mass() + 10 * p.width_at_pole();
        fdu = (umax - fumin) / N;
        fdm = (mmax - fmmin) / Nm;
        AddCoefficients(fmmin, mmax);
        TF1 f("f", [&](double *x, double *p) { return PDF(x[0]); }, fumin, umax, 0);
        for (int i = 0; i < N; ++i)
        {
            const double u = fumin + i * fdu;
            const double fmax = f.GetMaximum(u, u + fdu);
            fF.push_back(fmax);
            fu.push_back(u);
        }
        fu.push_back(umax);
        fConstDis = std::piecewise_constant_distribution<double>(fu.begin(), fu.end(), fF.begin());
        boost::math::quadrature::tanh_sinh<double> integrator;
        auto FMult = [&](double u)
        {
            return PDF(u);
        };
        fMultiplicity = volume * 0.25 * p.spin_degeneracy() * integrator.integrate(FMult, fumin, umax) * p.width_at_pole() / (constants::piSqr * constants::hbar3);
        double norm;
        if (fSF == SpectralFunction::BW)
            norm = (umax - fumin) / constants::pi;
        else
        {
            auto FNorm = [&](double u)
            {
                const double tanu = std::tan(u);
                const double m = p.mass() + 0.5 * p.width_at_pole() * tanu;
                return SpectralFunction(m) * (1.0 + tanu * tanu);
            };
            norm = integrator.integrate(FNorm, fumin, umax) * 0.5 * p.width_at_pole();
        }
        fMultiplicity /= norm;
    }
    fPoisson = std::poisson_distribution<int>(fMultiplicity);
}

std::tuple<double, double, double, double> ParticleSampler::operator()()
{
    int i;
    double m;
    // Sample particle mass
    if (fSF == SpectralFunction::None || fParticle.is_stable())
    {
        m = fParticle.mass();
        i = 0;
    }
    else
    {
        const double norm = 1.0 / fdu;
        double u = fConstDis(smash::random::engine);
        int j = static_cast<int>((u - fumin) * norm);
        if (fUni(smash::random::engine) * fF[j] < PDF(u))
            m = fParticle.mass() + 0.5 * fParticle.width_at_pole() * std::tan(u);
        else
        {
            do
            {
                u = fConstDis(smash::random::engine);
                j = static_cast<int>((u - fumin) * norm);
            } while (fUni(smash::random::engine) * fF[j] > PDF(u));
            m = fParticle.mass() + 0.5 * fParticle.width_at_pole() * std::tan(u);
        }
        if (m > fmmax)
            AddCoefficients(fmmax + fdm, m + fdm);
        i = (m - fmmin) / fdm;
    }
    // Given mass, sample particle momentum
    double p;
    do
        p = InverseMomentumCDF(fUni(smash::random::engine), i);
    while (fUni(smash::random::engine) * MomentumEnvelopeDF(p, i) > MomentumDistribution(p, m));
    const double phi = constants::twoPi * fUni(smash::random::engine);
    const double cost = 2.0 * fUni(smash::random::engine) - 1.0;
    const double sint = std::sqrt(1.0 - cost * cost);
    return {m, p * sint * std::cos(phi), p * sint * std::sin(phi), p * cost};
}

double ParticleSampler::MomentumEnvelopeDF(const double p, const int i) const
{
    if (p <= fX1[i])
        return fSlope[i] * p;
    if (p > fX2[i])
        return fMax[i] * std::exp(-fLambda[i] * (p - fX2[i]));
    return fMax[i];
}

double ParticleSampler::InverseMomentumCDF(const double u, const int i) const
{
    const double F1 = 0.5 * fSlope[i] * fX1[i] * fX1[i];
    const double F2 = F1 + fMax[i] * (fX2[i] - fX1[i]);
    const double CDFNorm = F2 + fMax[i] / fLambda[i];
    const double x = u * CDFNorm;
    if (x <= F1)
        return std::sqrt(2.0 * x / fSlope[i]);
    if (x > F2)
        return fX2[i] - std::log((CDFNorm - x) * fLambda[i] / fMax[i]) / fLambda[i];
    return (x - 0.5 * fSlope[i] * fX1[i] * fX1[i]) / fMax[i] + fX1[i];
}

double ParticleSampler::MomentumDistribution(const double p, const double m) const
{
    const double E = std::hypot(p, m);
    return p * p / (std::exp((E - fmu) / fT) + fa);
}

void ParticleSampler::AddCoefficients(const double mmin, const double mmax)
{
    for (double m = mmin; m <= mmax; m += fdm)
    {
        fmmax = m;
        auto FMax = [&](double x)
        {
            return 2.0 * std::hypot(x, m) * fT - x * x + fa * MomentumDistribution(x, m);
        };
        const auto &[xl1, xh1] = FindStartingPosition(0.0, 0.1, FMax);
        const double xMax = FindRoot(xl1, xh1, FMax);
        const double maxF = MomentumDistribution(xMax, m);
        fMax.push_back(maxF);

        auto FSlope = [&](double x)
        {
            return std::hypot(x, m) * fT - x * x + fa * MomentumDistribution(x, m);
        };
        const double xSlope = FindRoot(0.0, xMax, FSlope);
        const double slope = MomentumDistribution(xSlope, m) / xSlope;
        fX1.push_back(maxF / slope);
        fSlope.push_back(slope);

        auto FLog = [&](double x)
        {
            return std::log(MomentumDistribution(x, m)) - std::log(maxF) + 1.0;
        };
        const auto &[xl2, xh2] = FindStartingPosition(xMax, 0.1, FLog);
        const double xe = FindRoot(xl2, xh2, FLog);

        const double lambda = -(2.0 / xe - xe / (std::hypot(xe, m) * fT) + fa * MomentumDistribution(xe, m) / (xe * std::hypot(xe, m) * fT));
        const double x2 = xe - 1.0 / lambda;
        fX2.push_back(x2);
        fLambda.push_back(lambda);
    }
}

double ParticleSampler::SpectralFunction(const double m) const
{
    if (fSF == SpectralFunction::BW)
        return fParticle.spectral_function_simple(m);
    else if (fSF == SpectralFunction::RBW)
        return fParticle.spectral_function_const_width(m);
    else if (fSF == SpectralFunction::SMASH)
        return fParticle.spectral_function(m);
    return 0.0;
}

double ParticleSampler::PDF(const double u) const
{
    if (u < fumin || fSF == SpectralFunction::None)
        return 0.0;
    const double tanu = std::tan(u);
    const double m = fParticle.mass() + 0.5 * fParticle.width_at_pole() * tanu;
    const double n = Density(m);
    return SpectralFunction(m) * (1.0 + tanu * tanu) * n;
}

double ParticleSampler::Density(const double m) const
{
    const double mOverT = m / fT;
    const double exponent0 = std::exp(fmu / fT);
    double exponent = exponent0;
    double result = 0.0;
    double correction = 0.0;
    int k = 1;
    int ak = 1;
    double bk;
    do
    {
        try
        {
            bk = std::cyl_bessel_k(2, k * mOverT);
        }
        catch (...)
        {
            bk = 0.0;
        }
        correction = ak * bk * exponent / k;
        result += correction;
        exponent *= exponent0;
        ak *= -fa;
        ++k;
    } while (std::abs(correction / result) > 1e-6);
    return m * m * fT * result;
}
