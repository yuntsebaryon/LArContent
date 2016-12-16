/**
 *  @file   larpandoracontent/LArObjects/LArShowerPfo.h
 *
 *  @brief  Header file for the lar pfo class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_PFO_H
#define LAR_SHOWER_PFO_H 1

#include "Objects/ParticleFlowObject.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/ObjectFactory.h"

#include <string>

namespace lar_content
{

/**
 *  @brief  lar pfo parameters
 */
class LArShowerPfoParameters : public object_creation::ParticleFlowObject::Parameters
{
public:
    std::string     m_additionalProperty;       ///< The additional property string
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object
 */
class LArShowerPfo : public object_creation::ParticleFlowObject::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar pfo parameters
     */
    LArShowerPfo(const LArShowerPfoParameters &parameters);

    /**
     *  @brief  Get the additional property string
     * 
     *  @return the additional property string
     */
    const std::string &GetAdditionalProperty() const;

    /**
     *  @brief  Get the shower length from 3d shower fit
     * 
     *  @return the shower length from 3d shower fit
     */
    float GetShowerLength() const;
    
    /**
     *  @brief  Get the shower min layer position from 3d shower fit
     * 
     *  @return the shower min layer position from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerMinLayerPosition() const;

    /**
     *  @brief  Get the shower max layer position from 3d shower fit
     * 
     *  @return the shower max layer position from 3d shower fit
     */
    const pandora::CartesianVector &GetShowerMaxLayerPosition() const;

private:
    std::string                 m_additionalProperty;       ///< The additional property string

    float                       m_showerLength;             ///< Dummy property, shower length from 3d shower fit
    pandora::CartesianVector    m_showerMinLayerPosition;   ///< Dummy property, shower min layer position from 3d shower fit
    pandora::CartesianVector    m_showerMaxLayerPosition;   ///< Dummy property, shower max layer position from 3d shower fit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object factory responsible for pfo creation
 */
class LArShowerPfoFactory : public pandora::ObjectFactory<object_creation::ParticleFlowObject::Parameters, object_creation::ParticleFlowObject::Object>
{
public:
    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  fileReader the file reader, used to extract any additional parameters from file
     */
    pandora::StatusCode Read(Parameters &parameters, pandora::FileReader &fileReader) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const Object *const pObject, pandora::FileWriter &fileWriter) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const Parameters &parameters, const Object *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArShowerPfo::LArShowerPfo(const LArShowerPfoParameters &parameters) :
    object_creation::ParticleFlowObject::Object(parameters),
    m_additionalProperty(parameters.m_additionalProperty),
    m_showerLength(0),
    m_showerMinLayerPosition(0.f, 0.f, 0.f),
    m_showerMaxLayerPosition(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &LArShowerPfo::GetAdditionalProperty() const
{
    return m_additionalProperty;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArShowerPfoFactory::Parameters *LArShowerPfoFactory::NewParameters() const
{
    return (new LArShowerPfoParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArShowerPfoParameters &larPfoParameters(dynamic_cast<const LArShowerPfoParameters&>(parameters));
    pObject = new LArShowerPfo(larPfoParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Read(Parameters&, pandora::FileReader&) const
{
    // TODO: Provide this functionality if necessary

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArShowerPfoFactory::Write(const pandora::ParticleFlowObject*, pandora::FileWriter&) const
{
    // TODO: Provide this functionality if necessary

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_PFO_H
