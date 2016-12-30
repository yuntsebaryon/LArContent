/**
 *  @file   larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D shower building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArShowerPfo.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h"

#include "TPrincipal.h"

using namespace pandora;

namespace lar_content
{

ShowerParticleBuildingAlgorithm::ShowerParticleBuildingAlgorithm() :
    m_cosmicMode(false)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject*& pOutputPfo) const
{
    try
    {
        // Need an input vertex to provide a shower propagation direction
        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        // In cosmic mode, build showers from all daughter pfos, otherwise require that pfo is shower-like
        if (m_cosmicMode)
        {
            if (LArPfoHelper::IsFinalState(pInputPfo))
                return;
        }
        else
        {
            if (!LArPfoHelper::IsShower(pInputPfo))
                return;
        }

        // Build a new pfo
        LArShowerPfoFactory pfoFactory;
        LArShowerPfoParameters pfoParameters;
        pfoParameters.m_particleId = (LArPfoHelper::IsShower(pInputPfo) ? pInputPfo->GetParticleId() : E_MINUS);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = pInputPfo->GetMomentum();

        pfoParameters.m_showerLength = CartesianVector(0.f, 0.f, 0.f); // DUMMY
        pfoParameters.m_showerMinLayerPosition = CartesianVector(0.f, 0.f, 0.f); // DUMMY
        pfoParameters.m_showerMaxLayerPosition = CartesianVector(0.f, 0.f, 0.f); // DUMMY
        pfoParameters.m_showerCentroid = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_showerPositionSigmas = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_showerOpeningAngle = 0.f;
        pfoParameters.m_showerDirection = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_showerSecondaryVector = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_showerTertiaryVector = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_showerVertex = pInputVertex->GetPosition();

        ClusterList threeDClusterList;
        LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

        if (!threeDClusterList.empty())
        {
            const Cluster *const pThreeDCluster(threeDClusterList.front());

            const ClusterList listForVisualization(1, pThreeDCluster);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &listForVisualization, "ThreeDClusters", RED);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            CaloHitList threeDCaloHitList; // Might be able to get through LArPfoHelper instead
            pThreeDCluster->GetOrderedCaloHitList().FillCaloHitList(threeDCaloHitList);
            TPrincipal* principal = new TPrincipal( 3, "D" );

            for (const CaloHit *const pCaloHit3D : threeDCaloHitList)
            {
                double data[3];
                data[0] = pCaloHit3D->GetPositionVector().GetX();
                data[1] = pCaloHit3D->GetPositionVector().GetY();
                data[2] = pCaloHit3D->GetPositionVector().GetZ();
                principal->AddRow( data );
            }

            // Check if there are mean values
            if ( principal->GetMeanValues()->GetNrows() < 3 ) {
                // Need to discuss with John for the exception
                std::cerr << "Mean value issue!" << std::endl;
            }

            // Do the actual analysis
            principal->MakePrincipals();

            // Print out the result on
            principal->Print();

            try {
                pfoParameters.m_showerCentroid = CartesianVector( (*principal->GetMeanValues())[0], (*principal->GetMeanValues())[1], (*principal->GetMeanValues())[2] );
                pfoParameters.m_showerPositionSigmas = CartesianVector( (*principal->GetSigmas())[0], (*principal->GetSigmas())[1], (*principal->GetSigmas())[2] );
                pfoParameters.m_showerDirection = CartesianVector( (*principal->GetEigenVectors())[0][0], (*principal->GetEigenVectors())[1][0], (*principal->GetEigenVectors())[2][0] );
                pfoParameters.m_showerSecondaryVector = CartesianVector( (*principal->GetEigenVectors())[0][1], (*principal->GetEigenVectors())[1][1], (*principal->GetEigenVectors())[2][1] );
                pfoParameters.m_showerTertiaryVector = CartesianVector( (*principal->GetEigenVectors())[0][2], (*principal->GetEigenVectors())[1][2], (*principal->GetEigenVectors())[2][2] );
                const CartesianVector eigenvalues( (*principal->GetEigenValues())[0], (*principal->GetEigenValues())[1], (*principal->GetEigenValues())[2] );
                float norm = PCANormalization( pfoParameters.m_showerPositionSigmas.Get() );
                pfoParameters.m_showerEigenValues = eigenvalues * norm;
                CartesianVector sLength(0.f, 0.f, 0.f);
                ShowerLength( pfoParameters.m_showerEigenValues.Get(), sLength );
                pfoParameters.m_showerLength = sLength;
                pfoParameters.m_showerOpeningAngle = OpeningAngle( pfoParameters.m_showerDirection.Get(), pfoParameters.m_showerSecondaryVector.Get(), pfoParameters.m_showerEigenValues.Get() );

            } catch (const StatusCodeException &statusCodeException) {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }

/*
            const unsigned int layerHalfWindow(20); // TODO, read via xml
            const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

            try
            {
                const ThreeDSlidingFitResult threeDFitResult(pThreeDCluster, layerHalfWindow, layerPitch);

                const CartesianVector &globalMinLayerPosition(threeDFitResult.GetGlobalMinLayerPosition());
                const CartesianVector &globalMaxLayerPosition(threeDFitResult.GetGlobalMaxLayerPosition());

                pfoParameters.m_showerLength = ((globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude());
                pfoParameters.m_showerMinLayerPosition = globalMinLayerPosition;
                pfoParameters.m_showerMaxLayerPosition = globalMaxLayerPosition;

                std::cout << "ShowerLength " << pfoParameters.m_showerLength.Get() << std::endl;
            }
            catch (const StatusCodeException &)
            {
            } */
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pOutputPfo,
            pfoFactory));

        const LArShowerPfo *const pLArPfo = dynamic_cast<const LArShowerPfo*>(pOutputPfo);
        if (NULL == pLArPfo)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new vertex
        const Vertex *pOutputVertex(NULL);

        PandoraContentApi::Vertex::Parameters vtxParameters;
        vtxParameters.m_position = pInputVertex->GetPosition();
        vtxParameters.m_vertexLabel = pInputVertex->GetVertexLabel();
        vtxParameters.m_vertexType = pInputVertex->GetVertexType();

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vtxParameters, pOutputVertex));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pOutputPfo, pOutputVertex));
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
float ShowerParticleBuildingAlgorithm::PCANormalization( const pandora::CartesianVector Sigmas ) const
{
    return Sigmas.GetX() * Sigmas.GetX() + Sigmas.GetY() * Sigmas.GetY() + Sigmas.GetZ() * Sigmas.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------
void ShowerParticleBuildingAlgorithm::ShowerLength( const pandora::CartesianVector EigenValues, pandora::CartesianVector &sLength ) const
{
    sLength.SetValues( 6.*sqrt( EigenValues.GetX()), 6.*sqrt( EigenValues.GetY()), 6.*sqrt( EigenValues.GetZ()) );
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
float ShowerParticleBuildingAlgorithm::OpeningAngle( const pandora::CartesianVector principal, const pandora::CartesianVector secondary, const pandora::CartesianVector EigenValues ) const
{
    float cosTheta = principal.GetDotProduct( secondary ) / ( principal.GetMagnitude() * secondary.GetMagnitude() );
    float sinTheta = sqrt( 1. - cosTheta* cosTheta );
    float openAngle = 2.* atan( sqrt( EigenValues.GetY() ) * sinTheta / sqrt( EigenValues.GetX() ) );
    return openAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosmicMode", m_cosmicMode));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
