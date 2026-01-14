// C++
#include <iostream>
// Boost
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/beta.hpp>
// Project
#include "Constants.h"
#include "BlastWaveModus.h"
// SMASH
#include "smash/propagation.h"
#include "smash/configuration.h"
#include "smash/key.h"

BlastWaveModus::BlastWaveModus(smash::Configuration modusConfig, const smash::ExperimentParameters &param) : smash::ListModus()
{
    fUni = std::uniform_real_distribution<double>(0.0, 1.0);
    smash::Key<double> keyR({"Modi", "BlastWave", "R"}, {"1.0"});
    smash::Key<double> keyTau({"Modi", "BlastWave", "Tau"}, {"1.0"});
    smash::Key<double> keyEta({"Modi", "BlastWave", "Max_Eta"}, {"1.0"});
    smash::Key<double> keyV({"Modi", "BlastWave", "v"}, {"1.0"});
    smash::Key<double> keyN({"Modi", "BlastWave", "n"}, {"1.0"});
    smash::Key<double> keyT({"Modi", "BlastWave", "T"}, {"1.0"});
    smash::Key<double> keyMuB({"Modi", "BlastWave", "MuB"}, {"1.0"});
    smash::Key<double> keyMuS({"Modi", "BlastWave", "MuS"}, {"1.0"});
    smash::Key<double> keyMuQ({"Modi", "BlastWave", "MuQ"}, {"1.0"});
    smash::Key<double> keyMuPion({"Modi", "BlastWave", "MuPion"}, {"1.0"});
    smash::Key<double> keyMuKaon({"Modi", "BlastWave", "MuKaon"}, {"1.0"});
    smash::Key<std::string> keySF({"Modi", "BlastWave", "Spectral_Function"}, {"1.0"});
    fR = modusConfig.take(keyR);
    fTau = modusConfig.take(keyTau);
    fEtaMax = modusConfig.take(keyEta);
    fv = modusConfig.take(keyV);
    fn = modusConfig.take(keyN);
    fT = modusConfig.take(keyT);
    fMuB = modusConfig.take(keyMuB);
    fMuS = modusConfig.take(keyMuS);
    fMuQ = modusConfig.take(keyMuQ);
    fMuPion = modusConfig.take(keyMuPion);
    fMuKaon = modusConfig.take(keyMuKaon);
    std::string spectralFunction = modusConfig.take(keySF);

    if (spectralFunction == "BW")
        ParticleSampler::SetSpectralFunction(ParticleSampler::SpectralFunction::BW);
    else if (spectralFunction == "RBW")
        ParticleSampler::SetSpectralFunction(ParticleSampler::SpectralFunction::RBW);
    else if (spectralFunction == "SMASH")
        ParticleSampler::SetSpectralFunction(ParticleSampler::SpectralFunction::SMASH);
    else if (spectralFunction == "None")
        ParticleSampler::SetSpectralFunction(ParticleSampler::SpectralFunction::None);
    else
    {
        std::cout << "Unknown spectral function option! Automatically set to None!\n";
        ParticleSampler::SetSpectralFunction(ParticleSampler::SpectralFunction::None);
    }

    // Convert MeV to GeV
    fT *= 0.001;
    fMuB *= 0.001;
    fMuS *= 0.001;
    fMuQ *= 0.001;
    fMuPion *= 0.001;
    fMuKaon *= 0.001;
    // Print model summary
    smash::logg[smash::LogArea::Experiment::id].info() << "Blast-Wave Model Parameters:\n" 
    << "    T = " << fT << " GeV\n"
    << "    μᴮ = " << fMuB << " GeV\n"
    << "    μˢ = " << fMuS << " GeV\n"
    << "    μꟴ = " << fMuQ << " GeV\n"
    << "    μᴾⁱᵒⁿ = " << fMuPion << " GeV\n"
    << "    μᴷᵃᵒⁿ = " << fMuKaon << " GeV\n"
    << "    τ = " << fTau << " fm/c\n"
    << "    R = " << fR << " fm\n"
    << "    βₜ = " << fv << '\n'
    << "    n = " << fn << '\n'
    << "    ηᵐᵃˣ = " << fEtaMax << '\n'
    << "    Spectral function: " << spectralFunction << '\n';
    // Initialise particle samplers
    InitParticleInformation();
}

void BlastWaveModus::InitParticleInformation()
{
    smash::logg[smash::LogArea::Experiment::id].info("Initialising particle samplers"); 
    const double effectiveVolume = fEtaMax * constants::twoPi * fTau * fR * fR / (fn * std::pow(fv, 2.0 / fn)) * boost::math::beta(1.0 / fn, 0.5, fv * fv);
    for (const auto &p : smash::ParticleType::list_all())
    {
        // Skip non-hadrons
        if (!p.is_hadron())
            continue;

        // Skip particles with non-zero number of c and b quarks
        if (p.pdgcode().is_heavy_flavor())
            continue;

        const double mu = fMuB * p.baryon_number() + fMuS * p.strangeness() + fMuQ * p.charge() + fMuPion * p.is_pion() + fMuKaon * p.pdgcode().is_kaon();
        fParticleSamplers.emplace_back(p, fT, mu, effectiveVolume, 1000, 1000);
        smash::logg[smash::LogArea::Experiment::id].info() << "Initialised particle sampler for PDG code " << p.pdgcode().string() << ": <N> = " << fParticleSamplers.back().GetMeanMultiplicity(); 
    }
}

double BlastWaveModus::initial_conditions(smash::Particles *particles, const smash::ExperimentParameters &)
{
    const double maxWeight = std::exp(std::atanh(fv));
    for (auto &sampler : fParticleSamplers)
    {
        const int N = sampler.SampleMultiplicity();
        for (int i = 0; i < N; ++i)
        {
            double t, rx, ry, rz, m, E, px, py, pz, weight;
            do
            {
                const double eta = (2.0 * fUni(smash::random::engine) - 1.0) * fEtaMax;
                const double phi = fUni(smash::random::engine) * constants::twoPi;
                const double xr = fUni(smash::random::engine);
                const double coshEta = std::cosh(eta);
                const double sinhEta = std::sinh(eta);
                const double cosPhi = std::cos(phi);
                const double sinPhi = std::sin(phi);
                const double r = fR * std::sqrt(xr);
                t = fTau * coshEta;
                rx = r * cosPhi;
                ry = r * sinPhi;
                rz = fTau * sinhEta;
                const double coshRho = 1.0 / std::sqrt(1 - fv * fv * std::pow(xr, fn));
                const double sinhRho = std::sqrt(coshRho * coshRho - 1.0);
                const double u0 = coshRho * coshEta;
                const double ux = sinhRho * cosPhi;
                const double uy = sinhRho * sinPhi;
                const double uz = coshRho * sinhEta;
                const auto &[m0, pxl, pyl, pzl] = sampler();
                m = m0;
                const double El = std::sqrt(m * m + pxl * pxl + pyl * pyl + pzl * pzl);
                const double udotp = ux * pxl + uy * pyl + uz * pzl;
                E = u0 * El + udotp;
                const double factor = El + udotp / (1 + u0);
                px = pxl + ux * factor;
                py = pyl + uy * factor;
                pz = pzl + uz * factor;
                const double mt = std::hypot(m, px, py);
                const double y = 0.5 * std::log((E + pz) / (E - pz));
                weight = mt * std::cosh(y - eta) / El;
            } while (fUni(smash::random::engine) * maxWeight > weight);
            try_create_particle(*particles, sampler.GetPDGCode(), t, rx, ry, rz, m, E, px, py, pz);
        }
    }
    start_time_ = fTau;
    smash::backpropagate_straight_line(particles, start_time_);
    return start_time_;
}
