#ifndef _PARTICLE_SAMPLER_H_
#define _PARTICLE_SAMPLER_H_

// C++
#include <random>
#include <vector>
// SMASH
#include "smash/particletype.h"
#include "smash/random.h"


class ParticleSampler
{
public:
    enum class SpectralFunction
    {
        None,
        BW,
        RBW,
        SMASH
    };

    ParticleSampler(const smash::ParticleType &p, const double t, const double mu, const double volume, const int N = 1000, const int Nm = 1000);
    int SampleMultiplicity() { return fPoisson(smash::random::engine); }
    std::tuple<double, double, double, double> operator()();
    double MomentumEnvelopeDF(const double p, const int i) const;
    double InverseMomentumCDF(const double u, const int i) const;
    double GetMeanMultiplicity() const { return fMultiplicity; }
    bool IsStable() { return fParticle.is_stable(); }
    smash::PdgCode GetPDGCode() { return fParticle.pdgcode(); }
    const smash::ParticleType& GetParticleType() { return fParticle; }
    static void SetSpectralFunction(const SpectralFunction x) { fSF = x; }

private:
    const smash::ParticleType &fParticle;
    const double fT;
    const double fmu;
    double fdu;
    double fdm;
    double fumin;
    double fmmin;
    double fmmax;
    int fa;
    double fMultiplicity;

    std::vector<double> fu;
    std::vector<double> fF;
    std::vector<double> fMax;
    std::vector<double> fX1;
    std::vector<double> fX2;
    std::vector<double> fSlope;
    std::vector<double> fLambda;
    std::poisson_distribution<int> fPoisson;
    std::uniform_real_distribution<double> fUni;
    std::piecewise_constant_distribution<double> fConstDis;
    inline static SpectralFunction fSF = SpectralFunction::BW;

    void AddCoefficients(const double mmin, const double mmax);
    double MomentumDistribution(const double p, const double m) const;
    double SpectralFunction(const double m) const;
    double PDF(const double u) const;
    double Density(const double m) const;
    
    template <class F>
    std::pair<double, double> FindStartingPosition(const double x0, const double dx, const F& f, const int maxIt = 100)
    {
        double xl = x0;
        double xh = x0 + dx;
        int it = 0;
        while (f(xl) * f(xh) > 0.0 && it != maxIt)
        {
            xl += dx;
            xh += dx;
            ++it;
        }
        if (it == maxIt)
            throw std::runtime_error("Cannot find starting position for bisection! Reached iteration number limit!\n");
        return {xl, xh};
    }

    template <class F>
    double FindRoot(double xl, double xh, const F& f, const int maxIt = 100, const double tol = 1e-10)
    {
        double xm = 0.5 * (xl + xh);
        bool positiveL = (f(xl) > 0);
        int it = 0;
        do
        {
            const double fm = f(xm);
            if (positiveL && fm > 0.0)
                xl = xm;
            else
                xh = xm;
            xm = 0.5 * (xl + xh);
            ++it;
        } while (std::abs((xm - xl) / xm) > tol && it != maxIt);
        if (it == maxIt)
            std::cout << "Reached iteration number limit in bisection solver! Reached accuracy: " << std::abs((xm - xl) / xm) << ", f(x) = " << f(xm) << '\n';
        return xm;
    }
};

#endif