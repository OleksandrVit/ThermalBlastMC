#ifndef _Blast_Wave_Modus_H_
#define _Blast_Wave_Modus_H_

// C++
#include <vector>
// SMASH
#include "smash/listmodus.h"
// Project
#include "ParticleSampler.h"

class BlastWaveModus : public smash::ListModus
{
private:
    double fTau;
    double fR;
    double fEtaMax;
    double fT;
    double fMuB;
    double fMuS;
    double fMuQ;
    double fMuPion;
    double fMuKaon;
    double fn;
    double fv;
    std::vector<ParticleSampler> fParticleSamplers;
    std::uniform_real_distribution<double> fUni;
    void InitParticleInformation();

public:
    BlastWaveModus(smash::Configuration modusConfig, const smash::ExperimentParameters &param);
    // ~BlastWaveModus();
    double initial_conditions(smash::Particles *particles, const smash::ExperimentParameters &);
};

#endif // _FREEZE_OUT_MODEL_H_