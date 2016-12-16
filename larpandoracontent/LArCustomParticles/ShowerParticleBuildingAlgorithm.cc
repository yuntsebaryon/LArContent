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
#include "TRandom.h"

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


        // ATTN Dummy, test code
        const int n=10, m=1000;
        const int c = n / 5 + 1;

        TPrincipal* principal = new TPrincipal(n,"ND");

        // Use a pseudo-random number generator
        TRandom* random = new TRandom;

        // Make the m data-points
        // Make a variable to hold our data
        // Allocate memory for the data point
        Double_t* data = new Double_t[n];
        for (Int_t i = 0; i < m; i++) {

        // First we create the un-correlated, random variables, according
        // to one of three distributions
        for (Int_t j = 0; j < n - c; j++) {
           if (j % 3 == 0)
              data[j] = random->Gaus(5,1);
           else if (j % 3 == 1)
              data[j] = random->Poisson(8);
           else
              data[j] = random->Exp(2);
        }

        // Then we create the correlated variables
        for (Int_t j = 0 ; j < c; j++) {
           data[n - c + j] = 0;
           for (Int_t k = 0; k < n - c - j; k++)
              data[n - c + j] += data[k];
        }

        // Finally we're ready to add this datapoint to the PCA
        principal->AddRow(data);
        }

        // We delete the data after use, since TPrincipal got it by now.
        delete [] data;

        // Do the actual analysis
        principal->MakePrincipals();

        // Print out the result on
        principal->Print();

        // Build a new pfo
        LArShowerPfoFactory pfoFactory;
        LArShowerPfoParameters pfoParameters;
        pfoParameters.m_particleId = (LArPfoHelper::IsShower(pInputPfo) ? pInputPfo->GetParticleId() : E_MINUS);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = pInputPfo->GetMomentum();

        pfoParameters.m_showerLength = 0.f; // DUMMY
        pfoParameters.m_showerMinLayerPosition = CartesianVector(0.f, 0.f, 0.f); // DUMMY
        pfoParameters.m_showerMaxLayerPosition = CartesianVector(0.f, 0.f, 0.f); // DUMMY

        ClusterList threeDClusterList;
        LArPfoHelper::GetThreeDClusterList(pInputPfo, threeDClusterList);

        if (!threeDClusterList.empty())
        {
            const Cluster *const pThreeDCluster(threeDClusterList.front());

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
            }
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

StatusCode ShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosmicMode", m_cosmicMode));

    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
