/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/MissingTrackSegmentTool.cc
 * 
 *  @brief  Implementation of the missing track segment tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/LongTracksTool.h"
#include "LArThreeDReco/LArTrackMatching/MissingTrackSegmentTool.h"

using namespace pandora;

namespace lar
{

bool MissingTrackSegmentTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector; ClusterMergeMap clusterMergeMap;
    this->FindTracks(pAlgorithm, overlapTensor, protoParticleVector, clusterMergeMap);

    pAlgorithm->CreateThreeDParticles(protoParticleVector);
    const bool particlesMade(!protoParticleVector.empty());

    this->PerformClusterMerges(pAlgorithm, clusterMergeMap);
    const bool mergesMade(!clusterMergeMap.empty());

    return (particlesMade || mergesMade);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor,
    ProtoParticleVector &protoParticleVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(elementList, usedClusters, iteratorList);

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            if (!this->PassesParticleChecks(pAlgorithm, *(*iIter), overlapTensor, usedClusters, clusterMergeMap))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.insert((*iIter)->GetClusterU());
            protoParticle.m_clusterListV.insert((*iIter)->GetClusterV());
            protoParticle.m_clusterListW.insert((*iIter)->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const TransverseOverlapResult::XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());
        const float shortSpan(std::min(xOverlap.GetXSpanU(), std::min(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float longSpan1(std::max(xOverlap.GetXSpanU(), std::max(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float longSpan2(((xOverlap.GetXSpanU() > shortSpan) && (xOverlap.GetXSpanU() < longSpan1)) ? xOverlap.GetXSpanU() :
            ((xOverlap.GetXSpanV() > shortSpan) && (xOverlap.GetXSpanV() < longSpan1)) ? xOverlap.GetXSpanV() : xOverlap.GetXSpanW());

        if ((xOverlap.GetXOverlapSpan() < std::numeric_limits<float>::epsilon()) || (longSpan1 < std::numeric_limits<float>::epsilon()))
            continue;

        if ((shortSpan / xOverlap.GetXOverlapSpan()) < m_minXOverlapFraction)
            continue;

        if ((longSpan2 / longSpan1) < m_minXOverlapFraction)
            continue;

        iteratorList.push_back(eIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::PassesParticleChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element,
    const TensorType &overlapTensor, ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const
{
    const Particle particle(element);

    if ((particle.m_shortMinX - particle.m_longMinX < m_minParticleXEndSeparation) && (particle.m_longMaxX - particle.m_shortMaxX < m_minParticleXEndSeparation))
        return false;

    ClusterList candidateClusters;
    this->GetCandidateClusters(pAlgorithm, particle, overlapTensor, usedClusters, candidateClusters);

    if (candidateClusters.empty())
        return false;

    SlidingFitResultMap slidingFitResultMap;
    this->GetSlidingFitResultMap(pAlgorithm, candidateClusters, slidingFitResultMap);

    if (slidingFitResultMap.empty())
        return false;

    SegmentOverlapMap segmentOverlapMap;
    this->GetSegmentOverlapMap(pAlgorithm, particle, slidingFitResultMap, segmentOverlapMap);

    if (segmentOverlapMap.empty())
        return false;

    return this->MakeDecisions(particle, slidingFitResultMap, segmentOverlapMap, usedClusters, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetCandidateClusters(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle, const TensorType &overlapTensor,
    const ClusterList &usedClusters, ClusterList &candidateClusters) const
{
    const ClusterList &clusterList((TPC_VIEW_U == particle.m_shortHitType) ? pAlgorithm->GetInputClusterListU() :
        (TPC_VIEW_V == particle.m_shortHitType) ? pAlgorithm->GetInputClusterListV() : pAlgorithm->GetInputClusterListW());

    const TensorType::ClusterNavigationMap &navigationMap((TPC_VIEW_U == particle.m_shortHitType) ? overlapTensor.GetClusterNavigationMapUV() :
        (TPC_VIEW_V == particle.m_shortHitType) ? overlapTensor.GetClusterNavigationMapVW() : overlapTensor.GetClusterNavigationMapWU());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);

        if (pCluster == particle.m_pShortCluster)
            continue;

        if (navigationMap.count(pCluster) || usedClusters.count(pCluster))
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsInCandidateCluster)
            continue;

        candidateClusters.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetSlidingFitResultMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterList &candidateClusterList,
    SlidingFitResultMap &slidingFitResultMap) const
{
    for (ClusterList::const_iterator iter = candidateClusterList.begin(), iterEnd = candidateClusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);

        try
        {
            const TwoDSlidingFitResult &slidingFitResult(pAlgorithm->GetCachedSlidingFitResult(pCluster));
            slidingFitResultMap[pCluster] = slidingFitResult;
            continue;
        }
        catch (StatusCodeException &)
        {
        }

        try
        {
            TwoDSlidingFitResult slidingFitResult;
            LArClusterHelper::LArTwoDSlidingFit(pCluster, pAlgorithm->GetSlidingFitWindow(), slidingFitResult);
            slidingFitResultMap[pCluster] = slidingFitResult;
            continue;
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetSegmentOverlapMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle,
    const SlidingFitResultMap &slidingFitResultMap, SegmentOverlapMap &segmentOverlapMap) const
{
    const TwoDSlidingFitResult &fitResult1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster1));
    const TwoDSlidingFitResult &fitResult2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster2));
    const float nPoints1(std::fabs(static_cast<float>(fitResult1.GetMaxLayer() - fitResult1.GetMinLayer())));
    const float nPoints2(std::fabs(static_cast<float>(fitResult2.GetMaxLayer() - fitResult2.GetMinLayer())));
    const float xPitch((particle.m_longMaxX - particle.m_longMinX) * 2.f / (nPoints1 + nPoints2));

    for (float x = particle.m_longMinX; x < particle.m_longMaxX; x += xPitch)
    {
        if ((x > particle.m_shortMinX) && (x < particle.m_shortMaxX))
            continue;

        try
        {
            CartesianVector fitVector1(0.f, 0.f, 0.f), fitVector2(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitPosition(x, true, fitVector1);
            fitResult2.GetGlobalFitPosition(x, true, fitVector2);
            const float prediction(LArGeometryHelper::MergeTwoPositions(particle.m_hitType1, particle.m_hitType2, fitVector1.GetZ(), fitVector2.GetZ()));

            for (SlidingFitResultMap::const_iterator iter = slidingFitResultMap.begin(), iterEnd = slidingFitResultMap.end(); iter != iterEnd; ++iter)
            {
                try
                {
                    CartesianVector fitVector(0.f, 0.f, 0.f), fitDirection(0.f, 0.f, 0.f);
                    iter->second.GetGlobalFitPosition(x, true, fitVector);
                    iter->second.GetGlobalFitDirection(x, true, fitDirection);

                    const float delta((prediction - fitVector.GetZ()) * fitDirection.GetX());
                    const float pseudoChi2(delta * delta);

                    SegmentOverlap &segmentOverlap(segmentOverlapMap[iter->first]);
                    ++segmentOverlap.m_nSamplingPoints;
                    segmentOverlap.m_pseudoChi2Sum += pseudoChi2;

                    if (pseudoChi2 < m_pseudoChi2Cut)
                        ++segmentOverlap.m_nMatchedSamplingPoints;
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::MakeDecisions(const Particle &particle, const SlidingFitResultMap &slidingFitResultMap,
    const SegmentOverlapMap &segmentOverlapMap, ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const
{
    bool shouldMakeParticle(false);

    for (SegmentOverlapMap::const_iterator iter = segmentOverlapMap.begin(), iterEnd = segmentOverlapMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(iter->first);
        const SegmentOverlap &segmentOverlap(iter->second);

        if ((segmentOverlap.m_nSamplingPoints < m_makePfoMinSamplingPoints) || (segmentOverlap.m_nMatchedSamplingPoints < m_makePfoMinMatchedSamplingPoints))
            continue;

        if ((static_cast<float>(segmentOverlap.m_nMatchedSamplingPoints) / static_cast<float>(segmentOverlap.m_nSamplingPoints)) < m_makePfoMinMatchedFraction)
            continue;

        if (!this->AreDirectionsConsistent(particle.m_pShortCluster, pCluster))
            continue;

        shouldMakeParticle = true;
        usedClusters.insert(pCluster);

        if ((segmentOverlap.m_pseudoChi2Sum / static_cast<float>(segmentOverlap.m_nSamplingPoints)) > m_mergeMaxChi2PerSamplingPoint)
            continue;

        SlidingFitResultMap::const_iterator fitIter = slidingFitResultMap.find(pCluster);

        if (slidingFitResultMap.end() == fitIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        float minX(std::numeric_limits<float>::max()), maxX(-std::numeric_limits<float>::max());
        fitIter->second.GetMinAndMaxX(minX, maxX);

        if ((minX < particle.m_longMinX - m_mergeXContainmentTolerance) || (maxX > particle.m_longMaxX + m_mergeXContainmentTolerance))
            continue;

        clusterMergeMap[particle.m_pShortCluster].insert(pCluster);
    }

    return shouldMakeParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::AreDirectionsConsistent(Cluster *const pClusterA, Cluster *const pClusterB) const
{
    const LArPointingCluster pointingClusterA(pClusterA);
    const LArPointingCluster pointingClusterB(pClusterB);

    LArPointingCluster::Vertex vertexA, vertexB;
    LArPointingClusterHelper::GetClosestVertices(pointingClusterA, pointingClusterB, vertexA, vertexB);

    float transverseAB(std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
    float longitudinalAB(-std::numeric_limits<float>::max()), longitudinalBA(-std::numeric_limits<float>::max());

    LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);
    LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

    if (std::min(transverseAB, transverseBA) > m_makePfoMaxImpactParameter)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::PerformClusterMerges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterMergeMap &clusterMergeMap) const
{
    for (ClusterMergeMap::const_iterator pIter = clusterMergeMap.begin(), pIterEnd = clusterMergeMap.end(); pIter != pIterEnd; ++pIter)
    {
        if (pIter->second.empty())
            continue;

        const HitType parentHitType(LArThreeDHelper::GetClusterHitType(pIter->first));
        const std::string clusterListName((TPC_VIEW_U == parentHitType) ? pAlgorithm->GetClusterListNameU() : (TPC_VIEW_V == parentHitType) ? pAlgorithm->GetClusterListNameV() : pAlgorithm->GetClusterListNameW());

        for (ClusterList::const_iterator dIter = pIter->second.begin(), dIterEnd = pIter->second.end(); dIter != dIterEnd; ++dIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pIter->first, *dIter, clusterListName, clusterListName));
        }

        pAlgorithm->UpdateUponDeletion(pIter->first);
        pAlgorithm->UpdateForNewCluster(pIter->first);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MissingTrackSegmentTool::Particle::Particle(const TensorType::Element &element)
{
    const TransverseOverlapResult::XOverlap &xOverlap(element.GetOverlapResult().GetXOverlap());

    m_shortHitType = ((xOverlap.GetXSpanU() < xOverlap.GetXSpanV()) && (xOverlap.GetXSpanU() < xOverlap.GetXSpanW())) ? TPC_VIEW_U :
        ((xOverlap.GetXSpanV() < xOverlap.GetXSpanU()) && (xOverlap.GetXSpanV() < xOverlap.GetXSpanW())) ? TPC_VIEW_V :
        ((xOverlap.GetXSpanW() < xOverlap.GetXSpanU()) && (xOverlap.GetXSpanW() < xOverlap.GetXSpanV())) ? TPC_VIEW_W : CUSTOM;

    if (CUSTOM == m_shortHitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pShortCluster = (TPC_VIEW_U == m_shortHitType) ? element.GetClusterU() : (TPC_VIEW_V == m_shortHitType) ? element.GetClusterV() : element.GetClusterW();
    m_pCluster1 = (TPC_VIEW_U == m_shortHitType) ? element.GetClusterV() : (TPC_VIEW_V == m_shortHitType) ? element.GetClusterU() : element.GetClusterU();
    m_pCluster2 = (TPC_VIEW_U == m_shortHitType) ? element.GetClusterW() : (TPC_VIEW_V == m_shortHitType) ? element.GetClusterW() : element.GetClusterV();
    m_shortMinX = (TPC_VIEW_U == m_shortHitType) ? xOverlap.GetUMinX() : (TPC_VIEW_V == m_shortHitType) ? xOverlap.GetVMinX() : xOverlap.GetWMinX();
    m_shortMaxX = (TPC_VIEW_U == m_shortHitType) ? xOverlap.GetUMaxX() : (TPC_VIEW_V == m_shortHitType) ? xOverlap.GetVMaxX() : xOverlap.GetWMaxX();
    m_longMinX = (TPC_VIEW_U == m_shortHitType) ? std::max(xOverlap.GetVMinX(), xOverlap.GetWMinX()) : (TPC_VIEW_V == m_shortHitType) ? std::max(xOverlap.GetUMinX(), xOverlap.GetWMinX()) : std::max(xOverlap.GetUMinX(), xOverlap.GetWMinX());
    m_longMaxX = (TPC_VIEW_U == m_shortHitType) ? std::min(xOverlap.GetVMaxX(), xOverlap.GetWMaxX()) : (TPC_VIEW_V == m_shortHitType) ? std::min(xOverlap.GetUMaxX(), xOverlap.GetWMaxX()) : std::min(xOverlap.GetUMaxX(), xOverlap.GetWMaxX());

    m_hitType1 = LArThreeDHelper::GetClusterHitType(m_pCluster1);
    m_hitType2 = LArThreeDHelper::GetClusterHitType(m_pCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTrackSegmentTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minXOverlapFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    m_minParticleXEndSeparation = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinParticleXEndSeparation", m_minParticleXEndSeparation));

    m_minCaloHitsInCandidateCluster = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsInCandidateCluster", m_minCaloHitsInCandidateCluster));

    m_pseudoChi2Cut = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    m_makePfoMinSamplingPoints = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinSamplingPoints", m_makePfoMinSamplingPoints));

    m_makePfoMinMatchedSamplingPoints = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinMatchedSamplingPoints", m_makePfoMinMatchedSamplingPoints));

    m_makePfoMinMatchedFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinMatchedFraction", m_makePfoMinMatchedFraction));

    m_makePfoMaxImpactParameter = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMaxImpactParameter", m_makePfoMaxImpactParameter));

    m_mergeMaxChi2PerSamplingPoint = 0.25f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergeMaxChi2PerSamplingPoint", m_mergeMaxChi2PerSamplingPoint));

    m_mergeXContainmentTolerance = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergeXContainmentTolerance", m_mergeXContainmentTolerance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
