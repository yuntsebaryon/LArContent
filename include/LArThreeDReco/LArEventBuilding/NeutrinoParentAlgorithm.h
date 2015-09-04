/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoParentAlgorithm.h
 * 
 *  @brief  Header file for the neutrino parent algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
#define LAR_NEUTRINO_PARENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoParentAlgorithm class
 */
class NeutrinoParentAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Slice class
     */
    class Slice
    {
    public:
        pandora::CaloHitList    m_caloHitListU;                     ///< The u calo hit list
        pandora::CaloHitList    m_caloHitListV;                     ///< The v calo hit list
        pandora::CaloHitList    m_caloHitListW;                     ///< The w calo hit list
    };

    typedef std::vector<Slice> SliceList;

    typedef std::vector<pandora::HitType> HitTypeList;
    typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

    HitTypeList                 m_hitTypeList;                      ///< The hit type list
    HitTypeToNameMap            m_caloHitListNames;                 ///< The hit type to calo hit list name map
    HitTypeToNameMap            m_clusterListNames;                 ///< The hit type to cluster list name map

    std::string                 m_caloHitListNameU;                 ///< The name of the u input calo hit list
    std::string                 m_caloHitListNameV;                 ///< The name of the v input calo hit list
    std::string                 m_caloHitListNameW;                 ///< The name of the w input calo hit list
        
    std::string                 m_clusterListNameU;                 ///< The name of the u working cluster list
    std::string                 m_clusterListNameV;                 ///< The name of the v working cluster list
    std::string                 m_clusterListNameW;                 ///< The name of the w working cluster list
        
    std::string                 m_clusteringAlgorithm;              ///< The name of the two dimensional clustering algorithm
    std::string                 m_listDeletionAlgorithm;            ///< The name of the list deletion algorithm
    std::string                 m_listMovingAlgorithm;              ///< The name of the list moving algorithm

    pandora::StringVector       m_twoDAlgorithms;                   ///< The names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_threeDAlgorithms;                 ///< The names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_threeDHitAlgorithms;              ///< The names of the three dimensional hit creation algorithms
    pandora::StringVector       m_vertexAlgorithms;                 ///< The names of the vertex reconstruction algorithms
    pandora::StringVector       m_mopUpAlgorithms;                  ///< The names of the mop-up algorithms
    pandora::StringVector       m_neutrinoAlgorithms;               ///< The names of the neutrino building algorithms
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoParentAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoParentAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
