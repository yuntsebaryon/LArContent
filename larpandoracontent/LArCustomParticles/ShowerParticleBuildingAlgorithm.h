/**
 *  @file   larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h
 *
 *  @brief  Header file for the neutrino event creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
#define LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArShowerPfo.h"

#include "larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerParticleBuildingAlgorithm class
 */
class ShowerParticleBuildingAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerParticleBuildingAlgorithm();

    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject*& pOutputPfo) const;

    float PCANormalization( const pandora::CartesianVector Sigmas ) const;

    void ShowerLength( const pandora::CartesianVector EigenValues, pandora::CartesianVector &ShowerLength ) const;

    float OpeningAngle( const pandora::CartesianVector principal, const pandora::CartesianVector secondary, const pandora::CartesianVector EigenValues ) const;

    bool            m_cosmicMode;             ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ShowerParticleBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ShowerParticleBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
