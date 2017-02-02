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

#include <Eigen/Dense>
// #include "TPrincipal.h"

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
            CartesianVector Centroid( 0.f, 0.f, 0.f ), EigenValues( 0.f, 0.f, 0.f );
            EigenVectors EigenVecs;
            RunPCA( threeDCaloHitList, Centroid, EigenValues, EigenVecs );

            try {
                pfoParameters.m_showerCentroid = Centroid;
                // pfoParameters.m_showerPositionSigmas;
                pfoParameters.m_showerDirection = EigenVecs[0];
                pfoParameters.m_showerSecondaryVector = EigenVecs[1];
                pfoParameters.m_showerTertiaryVector = EigenVecs[2];
                pfoParameters.m_showerEigenValues = EigenValues;
                pfoParameters.m_showerLength = ShowerLength( pfoParameters.m_showerEigenValues.Get() );
                pfoParameters.m_showerOpeningAngle = OpeningAngle( pfoParameters.m_showerDirection.Get(), pfoParameters.m_showerSecondaryVector.Get(), pfoParameters.m_showerEigenValues.Get() );

            } catch ( const StatusCodeException &statusCodeException ) {
                if ( STATUS_CODE_FAILURE == statusCodeException.GetStatusCode() )
                    throw statusCodeException;
            }


            const unsigned int layerHalfWindow(20); // TODO, read via xml
            const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

            try
            {
                const ThreeDSlidingFitResult threeDFitResult(pThreeDCluster, layerHalfWindow, layerPitch);

                const CartesianVector &globalMinLayerPosition(threeDFitResult.GetGlobalMinLayerPosition());
                const CartesianVector &globalMaxLayerPosition(threeDFitResult.GetGlobalMaxLayerPosition());

                pfoParameters.m_showerMinLayerPosition = globalMinLayerPosition;
                pfoParameters.m_showerMaxLayerPosition = globalMaxLayerPosition;
            }
            catch ( const StatusCodeException &statusCodeException )
            {
                if ( STATUS_CODE_FAILURE == statusCodeException.GetStatusCode() )
                    throw statusCodeException;
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
void ShowerParticleBuildingAlgorithm::RunPCA( const pandora::CaloHitList& threeDCaloHitList, pandora::CartesianVector& Centroid, pandora::CartesianVector& EigenValues, EigenVectors& EigenVecs ) const
{
    // The steps are:
    // 1) do a mean normalization of the input vec points
    // 2) compute the covariance matrix
    // 3) run the SVD
    // 4) extract the eigen vectors and values
    // see what happens

    // Run through the CaloHitList and get the mean position of all the hits
    double meanPos[] = { 0., 0., 0. };
    int    numthreeDHitsInt(0);

    for (const CaloHit *const pCaloHit3D : threeDCaloHitList)
    {
        meanPos[0] += pCaloHit3D->GetPositionVector().GetX();
        meanPos[1] += pCaloHit3D->GetPositionVector().GetY();
        meanPos[2] += pCaloHit3D->GetPositionVector().GetZ();
        numthreeDHitsInt++;
    }

    if ( numthreeDHitsInt == 0 ) {
        std::cout << "No three dimensional hit found!" << std::endl;
        throw StatusCodeException( STATUS_CODE_NOT_FOUND );
    }

    double numthreeDHits = double( numthreeDHitsInt );
    meanPos[0] /= numthreeDHits;
    meanPos[1] /= numthreeDHits;
    meanPos[2] /= numthreeDHits;
    Centroid = CartesianVector( meanPos[0], meanPos[1], meanPos[2] );

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.0);
    double yi2(0.0);
    double yizi(0.0);
    double zi2(0.);
    double weightSum(0.);

    for ( const CaloHit *const pCaloHit3D : threeDCaloHitList )
    {
        double weight(1.);
        double x = ( pCaloHit3D->GetPositionVector().GetX() - meanPos[0] ) * weight;
        double y = ( pCaloHit3D->GetPositionVector().GetY() - meanPos[1] ) * weight;
        double z = ( pCaloHit3D->GetPositionVector().GetZ() - meanPos[2] ) * weight;

        weightSum += weight*weight;

        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
    }

    // Using Eigen package
    Eigen::Matrix3f sig;

    sig <<  xi2, xiyi, xizi,
           xiyi,  yi2, yizi,
           xizi, yizi,  zi2;

    if ( weightSum == 0. ) {
        std::cout << "The total weight of three dimensional hits is 0!" << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }

    sig *= 1./weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

    if ( eigenMat.info() == Eigen::ComputationInfo::Success )
    {
        using eigenValColPair = std::pair<float,size_t>;
        std::vector<eigenValColPair> eigenValColVec;

        eigenValColVec.push_back(eigenValColPair(eigenMat.eigenvalues()(0),0));
        eigenValColVec.push_back(eigenValColPair(eigenMat.eigenvalues()(1),1));
        eigenValColVec.push_back(eigenValColPair(eigenMat.eigenvalues()(2),2));

        std::sort( eigenValColVec.begin(), eigenValColVec.end(), [](const eigenValColPair& left, const eigenValColPair& right){return left.first > right.first;} );

        // Now copy output
        // Get the eigen values
        EigenValues = CartesianVector( eigenValColVec[0].first, eigenValColVec[1].first, eigenValColVec[2].first );

        // Grab the principle axes
        Eigen::Matrix3f eigenVecs(eigenMat.eigenvectors());

        for ( const auto& pair : eigenValColVec )
        {
            CartesianVector tempVec = CartesianVector( eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second) );
            EigenVecs.push_back(tempVec);
        }
    }
    else
    {
        std::cout << "PCA decompose failure, number of three D hits = " << numthreeDHits << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
pandora::CartesianVector ShowerParticleBuildingAlgorithm::ShowerLength( const pandora::CartesianVector& EigenValues ) const
{
    double sl[] = { 0., 0., 0. };
    if ( EigenValues.GetX() > 0. ) sl[0] = 6.*std::sqrt( EigenValues.GetX() );
    else {
        std::cout << "The principal eigenvalue is equal to or less than 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }
    if ( EigenValues.GetY() > 0. ) sl[1] = 6.*std::sqrt( EigenValues.GetY() );
    if ( EigenValues.GetZ() > 0. ) sl[2] = 6.*std::sqrt( EigenValues.GetZ() );
    CartesianVector sLength( sl[0], sl[1], sl[2] );
    return sLength;
}

//------------------------------------------------------------------------------------------------------------------------------------------
float ShowerParticleBuildingAlgorithm::OpeningAngle( const pandora::CartesianVector& principal, const pandora::CartesianVector& secondary, const pandora::CartesianVector& EigenValues ) const
{
    if ( principal.GetMagnitude() == 0. ) {
        std::cout << "The principal eigenvector is 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    } else if ( secondary.GetMagnitude() == 0. ) return 0.;

    float cosTheta = principal.GetDotProduct( secondary ) / ( principal.GetMagnitude() * secondary.GetMagnitude() );
    if ( cosTheta > 1. ) {
        std::cout << "cos(theta) between the principal and secondary eigenvectors is greater than 1." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    }

    float sinTheta = std::sqrt( 1. - cosTheta* cosTheta );
    if ( EigenValues.GetX() <= 0. ) {
        std::cout << "The principal eigenvalue is equal to or less than 0." << std::endl;
        throw StatusCodeException( STATUS_CODE_INVALID_PARAMETER );
    } else if ( EigenValues.GetY() < 0. ) return 0.;

    float openAngle = 2.* std::atan( std::sqrt( EigenValues.GetY() ) * sinTheta / std::sqrt( EigenValues.GetX() ) );
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
