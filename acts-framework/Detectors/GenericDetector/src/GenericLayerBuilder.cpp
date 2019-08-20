// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/GenericDetector/GenericLayerBuilder.hpp"
#include <iostream>
#include "ACTFW/GenericDetector/GenericDetectorElement.hpp"
#include "Acts/Detector/DetectorElementBase.hpp"
#include "Acts/Layers/ProtoLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Plugins/Digitization/CartesianSegmentation.hpp"
#include "Acts/Plugins/Digitization/DigitizationModule.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Tools/LayerCreator.hpp"
#include "Acts/Utilities/ApproachDescriptor.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

FW::Generic::GenericLayerBuilder::GenericLayerBuilder(
    const FW::Generic::GenericLayerBuilder::Config& glbConfig,
    std::unique_ptr<const Acts::Logger>             log)
  : Acts::ILayerBuilder()
  , m_nLayers()
  , m_cLayers()
  , m_pLayers()
  , m_logger(std::move(log))
{
  /// @todo a configuraiton check should be done here
  setConfiguration(glbConfig);
  // Tool needs to be initialized
  constructLayers();
}

void
FW::Generic::GenericLayerBuilder::setConfiguration(
    const FW::Generic::GenericLayerBuilder::Config& glbConfig)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = glbConfig;
}

void
FW::Generic::GenericLayerBuilder::setLogger(
    std::unique_ptr<const Acts::Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

void
FW::Generic::GenericLayerBuilder::constructLayers()
{
  size_t imodule = 0;
  // ----------------------- central layers -------------------------
  // the central layers
  size_t numcLayers = m_cfg.centralLayerRadii.size();
  if (numcLayers) {
    ACTS_DEBUG("Configured to build " << numcLayers
                                      << " active central layers.");
    m_cLayers.reserve(numcLayers);
    // loop through
    for (size_t icl = 0; icl < numcLayers; ++icl) {
      // layer R/Z
      double layerR = m_cfg.centralLayerRadii.at(icl);
      // some screen output
      ACTS_DEBUG("Build layer " << icl << " with target radius = " << layerR);

      // prepare the Surface vector
      std::vector<std::shared_ptr<const Acts::Surface>> sVector;
      // assign the current envelope
      double layerEnvelopeCoverZ = m_cfg.centralLayerEnvelopes.size()
          ? m_cfg.centralLayerEnvelopes.at(icl).second
          : 0.;
      // module size & tilt
      double modulePhiTilt   = m_cfg.centralModuleTiltPhi.at(icl);
      double moduleHalfX     = m_cfg.centralModuleHalfX.at(icl);
      double moduleHalfY     = m_cfg.centralModuleHalfY.at(icl);
      double moduleThickness = m_cfg.centralModuleThickness.at(icl);
      // create the shared module
      std::shared_ptr<const Acts::PlanarBounds> moduleBounds(
          new Acts::RectangleBounds(moduleHalfX, moduleHalfY));
      // Identifier @todo unique Identifier - use a GenericDetector identifier
      size_t nCetralModules = m_cfg.centralModuleBinningSchema.at(icl).first
          * m_cfg.centralModuleBinningSchema.at(icl).second;

      ACTS_DEBUG("- number of modules "
                 << nCetralModules
                 << " ( from "
                 << m_cfg.centralModuleBinningSchema.at(icl).first
                 << " x "
                 << m_cfg.centralModuleBinningSchema.at(icl).second
                 << " )");

      sVector.reserve(nCetralModules);

      // prepartion :
      // create digitizaiton module
      std::shared_ptr<const Acts::DigitizationModule> moduleDigitizationPtr
          = nullptr;
      if (m_cfg.centralModuleReadoutBinsX.size()) {
        // create the CartesianSegmentation
        std::shared_ptr<const Acts::Segmentation> moduleSegmentation
            = std::make_shared<const Acts::CartesianSegmentation>(
                moduleBounds,
                m_cfg.centralModuleReadoutBinsX.at(icl),
                m_cfg.centralModuleReadoutBinsY.at(icl));
        // now create the digitzation module
        moduleDigitizationPtr
            = std::make_shared<const Acts::DigitizationModule>(
                moduleSegmentation,
                m_cfg.centralModuleThickness.at(icl),
                m_cfg.centralModuleReadoutSide.at(icl),
                m_cfg.centralModuleLorentzAngle.at(icl));
      }

      // prepartation :
      // create the Module material from input
      std::shared_ptr<const Acts::SurfaceMaterial> moduleMaterialPtr = nullptr;
      if (m_cfg.centralModuleMaterial.size()) {
        // get the sensor material from configuration
        Acts::Material moduleMaterial = m_cfg.centralModuleMaterial.at(icl);
        Acts::MaterialProperties moduleMaterialProperties(moduleMaterial,
                                                          moduleThickness);
        // create a new surface material
        moduleMaterialPtr = std::shared_ptr<const Acts::SurfaceMaterial>(
            new Acts::HomogeneousSurfaceMaterial(moduleMaterialProperties));
      }

      // confirm
      if (m_cfg.centralModulePositions.at(icl).size() != nCetralModules) {
        ACTS_WARNING("Mismatching module numbers, configuration error!");
        ACTS_WARNING("- Binning schema suggests : " << nCetralModules);
        ACTS_WARNING("- Positions provided are  : "
                     << m_cfg.centralModulePositions.at(icl).size());
      }
      // loop over the position, create the modules
      for (auto& moduleCenter : m_cfg.centralModulePositions.at(icl)) {
        // create the association transform
        double modulePhi = phi(moduleCenter);
        // the local z axis is the normal vector
        Acts::Vector3D moduleLocalZ(
            cos(modulePhi + modulePhiTilt), sin(modulePhi + modulePhiTilt), 0.);
        // the local y axis is the global z axis
        Acts::Vector3D moduleLocalY(0., 0., 1);
        // the local x axis the normal to local y,z
        Acts::Vector3D moduleLocalX(-sin(modulePhi + modulePhiTilt),
                                    cos(modulePhi + modulePhiTilt),
                                    0.);
        // create the RotationMatrix
        Acts::RotationMatrix3D moduleRotation;
        moduleRotation.col(0) = moduleLocalX;
        moduleRotation.col(1) = moduleLocalY;
        moduleRotation.col(2) = moduleLocalZ;
        // get the moduleTransform
        std::shared_ptr<Acts::Transform3D> mutableModuleTransform
            = std::make_shared<Acts::Transform3D>(
                Acts::Translation3D(moduleCenter) * moduleRotation);
        // stereo angle if necessary
        if (m_cfg.centralModuleFrontsideStereo.size()
            && m_cfg.centralModuleFrontsideStereo.at(icl) != 0.) {
          // twist by the stereo angle
          double stereo = m_cfg.centralModuleFrontsideStereo.at(icl);
          (*mutableModuleTransform.get())
              *= Acts::AngleAxis3D(-stereo, Acts::Vector3D::UnitZ());
        }
        // count the modules
        ++imodule;
        Identifier moduleIdentifier
            = Identifier(Identifier::identifier_type(imodule));
        // Finalize the transform
        auto moduleTransform = std::const_pointer_cast<const Acts::Transform3D>(
            mutableModuleTransform);
        // create the module
        Acts::DetectorElementBase* module
            = new FW::Generic::GenericDetectorElement(moduleIdentifier,
                                                      moduleTransform,
                                                      moduleBounds,
                                                      moduleThickness,
                                                      moduleMaterialPtr,
                                                      moduleDigitizationPtr);
        // register the surface
        sVector.push_back(module->surface().getSharedPtr());
        // store the module
        // @todo detector store facility
        m_centralModule.push_back(module);
        // IF double modules exist
        // and the backside one (if configured to do so)
        if (m_cfg.centralModuleBacksideGap.size()) {
          // ncrease the counter @todo switch to identifier service
          ++imodule;
          // create the module identifier
          moduleIdentifier = Identifier(Identifier::identifier_type(imodule));
          moduleCenter     = moduleCenter
              + m_cfg.centralModuleBacksideGap.at(icl) * moduleLocalZ;
          mutableModuleTransform = std::make_shared<Acts::Transform3D>(
              Acts::Translation3D(moduleCenter) * moduleRotation);
          // apply the stereo
          if (m_cfg.centralModuleBacksideStereo.size()) {
            // twist by the stereo angle
            double stereoBackSide = m_cfg.centralModuleBacksideStereo.at(icl);
            (*mutableModuleTransform.get())
                *= Acts::AngleAxis3D(-stereoBackSide, Acts::Vector3D::UnitZ());
          }
          // Finalize the transform
          moduleTransform = std::const_pointer_cast<const Acts::Transform3D>(
              mutableModuleTransform);
          // everything is set for the next module
          Acts::DetectorElementBase* bsmodule
              = new FW::Generic::GenericDetectorElement(moduleIdentifier,
                                                        moduleTransform,
                                                        moduleBounds,
                                                        moduleThickness,
                                                        moduleMaterialPtr,
                                                        moduleDigitizationPtr);
          // register the backside as bin member
          std::vector<const Acts::DetectorElementBase*> bsbinmember = {module};
          std::vector<const Acts::DetectorElementBase*> binmember = {bsmodule};
          bsmodule->registerBinmembers(bsbinmember);
          module->registerBinmembers(binmember);
          // memory management - we need a detector store to hold them
          // somewhere @todo detector store facility
          m_centralModule.push_back(bsmodule);
        }
      }

      size_t phiBins = m_cfg.centralModuleBinningSchema.at(icl).first;
      phiBins *= m_cfg.centralLayerBinMultipliers.first;
      size_t zBins = m_cfg.centralModuleBinningSchema.at(icl).second;
      zBins *= m_cfg.centralLayerBinMultipliers.second;
      // create the surface array - it will also fill the accesible binmember
      // chache if avalable
      Acts::ProtoLayer pl(sVector);
      pl.envR = {m_cfg.approachSurfaceEnvelope, m_cfg.approachSurfaceEnvelope};
      pl.envZ = {layerEnvelopeCoverZ, layerEnvelopeCoverZ};
      Acts::MutableLayerPtr cLayer
          = m_cfg.layerCreator->cylinderLayer(sVector, phiBins, zBins, pl);
      // the layer is built le't see if it needs material
      if (m_cfg.centralLayerMaterialProperties.size()) {
        // get the material from configuration
        Acts::MaterialProperties layerMaterialProperties
            = m_cfg.centralLayerMaterialProperties.at(icl);
        std::shared_ptr<const Acts::SurfaceMaterial> layerMaterialPtr(
            new Acts::HomogeneousSurfaceMaterial(layerMaterialProperties));
        // central material
        if (m_cfg.centralLayerMaterialConcentration.at(icl) == 0.) {
          // the layer surface is the material surface
          cLayer->surfaceRepresentation().setAssociatedMaterial(
              layerMaterialPtr);
          ACTS_VERBOSE("- and material at central layer surface.");
        } else {
          // approach surface material
          // get the approach descriptor - at this stage we know that the
          // approachDescriptor exists
          auto approachSurfaces
              = cLayer->approachDescriptor()->containedSurfaces();
          if (m_cfg.centralLayerMaterialConcentration.at(icl) > 0) {
            auto mutableOuterSurface
                = const_cast<Acts::Surface*>(approachSurfaces.at(1));
            mutableOuterSurface->setAssociatedMaterial(layerMaterialPtr);
            ACTS_VERBOSE("- and material at outer approach surface");
          } else {
            auto mutableInnerSurface
                = const_cast<Acts::Surface*>(approachSurfaces.at(0));
            mutableInnerSurface->setAssociatedMaterial(layerMaterialPtr);
            ACTS_VERBOSE("- and material at inner approach surface");
          }
        }
      }
      // push it into the layer vector
      m_cLayers.push_back(cLayer);
    }
  }

  // -------------------------------- endcap type layers
  // pos/neg layers
  size_t numpnLayers = m_cfg.posnegLayerPositionsZ.size();
  if (numpnLayers) {
    ACTS_DEBUG("Configured to build 2 * "
               << numpnLayers
               << " passive positive/negative side layers.");
    m_pLayers.reserve(numpnLayers);
    m_nLayers.reserve(numpnLayers);

    /// this is the loop over th elayer positions
    for (size_t ipnl = 0; ipnl < numpnLayers; ++ipnl) {
      // some screen output
      ACTS_VERBOSE(
          "- building layers " << ipnl << " and " << numpnLayers + ipnl
                               << " at +/- z = "
                               << m_cfg.posnegLayerPositionsZ.at(ipnl));
      /// some preparation work
      // define the layer envelope
      double layerEnvelopeR = m_cfg.posnegLayerEnvelopeR.at(ipnl);
      // prepare for the r binning
      std::vector<std::shared_ptr<const Acts::Surface>> nsVector;
      std::vector<std::shared_ptr<const Acts::Surface>> psVector;
      // now fill the vectors
      size_t ipnR = 0;
      for (auto& discModulePositions : m_cfg.posnegModulePositions.at(ipnl)) {
        ACTS_VERBOSE("- building ring " << ipnR << " for this pair.");
        // now prepare all the shared stuff
        // (0) module specifications
        double moduleThickness = m_cfg.posnegModuleThickness.at(ipnl).at(ipnR);
        double moduleMinHalfX  = m_cfg.posnegModuleMinHalfX.at(ipnl).at(ipnR);
        double moduleMaxHalfX  = 0.;
        if (m_cfg.posnegModuleMaxHalfX.size() > ipnl
            && m_cfg.posnegModuleMaxHalfX.at(ipnl).size() > ipnR) {
          moduleMaxHalfX = m_cfg.posnegModuleMaxHalfX.at(ipnl).at(ipnR);
        }
        double moduleHalfY = m_cfg.posnegModuleHalfY.at(ipnl).at(ipnR);
        // (1) module bounds
        // create the bounds
        Acts::PlanarBounds* pBounds = nullptr;
        if (moduleMaxHalfX != 0. && moduleMinHalfX != moduleMaxHalfX)
          pBounds = new Acts::TrapezoidBounds(
              moduleMinHalfX, moduleMaxHalfX, moduleHalfY);
        else
          pBounds = new Acts::RectangleBounds(moduleMinHalfX, moduleHalfY);
        // now create the shared bounds from it
        std::shared_ptr<const Acts::PlanarBounds> moduleBounds(pBounds);
        // (2) create digitizaiton module
        std::shared_ptr<const Acts::DigitizationModule> moduleDigitizationPtr
            = nullptr;
        if (m_cfg.posnegModuleReadoutBinsX.size()) {
          // create the CartesianSegmentation
          std::shared_ptr<const Acts::Segmentation> moduleSegmentation
              = std::make_shared<const Acts::CartesianSegmentation>(
                  moduleBounds,
                  m_cfg.posnegModuleReadoutBinsX.at(ipnl).at(ipnR),
                  m_cfg.posnegModuleReadoutBinsY.at(ipnl).at(ipnR));
          // now create the digitzation module
          moduleDigitizationPtr
              = std::make_shared<const Acts::DigitizationModule>(
                  moduleSegmentation,
                  moduleThickness,
                  m_cfg.posnegModuleReadoutSide.at(ipnl).at(ipnR),
                  m_cfg.posnegModuleLorentzAngle.at(ipnl).at(ipnR));
        }
        // (3) module material
        // create the Module material from input
        std::shared_ptr<const Acts::SurfaceMaterial> moduleMaterialPtr
            = nullptr;
        if (m_cfg.posnegModuleMaterial.size()) {
          Acts::MaterialProperties moduleMaterialProperties(
              m_cfg.posnegModuleMaterial.at(ipnl).at(ipnR), moduleThickness);
          // and create the shared pointer
          moduleMaterialPtr = std::shared_ptr<const Acts::SurfaceMaterial>(
              new Acts::HomogeneousSurfaceMaterial(moduleMaterialProperties));
        }

        // low loop over the phi positions and build the stuff
        for (auto& ringModulePosition : discModulePositions) {
          // the module transform from the position
          double modulePhi = phi(ringModulePosition);
          // the center position of the modules
          Acts::Vector3D pModuleCenter(ringModulePosition);
          // take the mirrored position wrt x/y
          Acts::Vector3D nModuleCenter(
              pModuleCenter.x(), pModuleCenter.y(), -pModuleCenter.z());
          // the rotation matrix of the module
          Acts::Vector3D moduleLocalY(cos(modulePhi), sin(modulePhi), 0.);
          // take different axis to have the same readout direction
          Acts::Vector3D pModuleLocalZ(0., 0., 1.);
          // take different axis to have the same readout direction
          Acts::Vector3D nModuleLocalZ(0., 0., -1.);
          Acts::Vector3D nModuleLocalX = moduleLocalY.cross(nModuleLocalZ);
          Acts::Vector3D pModuleLocalX = moduleLocalY.cross(pModuleLocalZ);
          // local rotation matrices
          // create the RotationMatrix - negative side
          Acts::RotationMatrix3D nModuleRotation;
          nModuleRotation.col(0) = nModuleLocalX;
          nModuleRotation.col(1) = moduleLocalY;
          nModuleRotation.col(2) = nModuleLocalZ;
          // create the RotationMatrix - positive side
          Acts::RotationMatrix3D pModuleRotation;
          pModuleRotation.col(0) = pModuleLocalX;
          pModuleRotation.col(1) = moduleLocalY;
          pModuleRotation.col(2) = pModuleLocalZ;
          // the transforms for the two modules
          std::shared_ptr<const Acts::Transform3D> nModuleTransform
              = std::make_shared<const Acts::Transform3D>(
                  Acts::Translation3D(nModuleCenter) * nModuleRotation);
          std::shared_ptr<const Acts::Transform3D> pModuleTransform
              = std::make_shared<const Acts::Transform3D>(
                  Acts::Translation3D(pModuleCenter) * pModuleRotation);
          // create the modules identifier @todo Idenfier service
          Identifier nModuleIdentifier
              = Identifier(Identifier::identifier_type(2 * imodule));
          Identifier pModuleIdentifier
              = Identifier(Identifier::identifier_type(2 * imodule + 1));
          // create the module
          FW::Generic::GenericDetectorElement* nmodule
              = new FW::Generic::GenericDetectorElement(nModuleIdentifier,
                                                        nModuleTransform,
                                                        moduleBounds,
                                                        moduleThickness,
                                                        moduleMaterialPtr);
          FW::Generic::GenericDetectorElement* pmodule
              = new FW::Generic::GenericDetectorElement(pModuleIdentifier,
                                                        pModuleTransform,
                                                        moduleBounds,
                                                        moduleThickness,
                                                        moduleMaterialPtr);
          // memory management - we need a detector store to hold them somewhere
          // @todo add detector store facility
          m_posnegModule.push_back(nmodule);
          m_posnegModule.push_back(pmodule);
          // now deal with the potential backside
          if (m_cfg.posnegModuleBacksideGap.size()) {
            // ncrease the counter @todo switch to identifier service
            nModuleIdentifier
                = Identifier(Identifier::identifier_type(++imodule));
            pModuleIdentifier
                = Identifier(Identifier::identifier_type(++imodule));
            // the new centers
            nModuleCenter = nModuleCenter
                + m_cfg.posnegModuleBacksideGap.at(ipnl).at(ipnR)
                    * nModuleLocalZ;
            pModuleCenter = pModuleCenter
                + m_cfg.posnegModuleBacksideGap.at(ipnl).at(ipnR)
                    * pModuleLocalZ;
            // the new transforms
            auto mutableNModuleTransform = std::make_shared<Acts::Transform3D>(
                Acts::Translation3D(nModuleCenter) * nModuleRotation);
            auto mutablePModuleTransform = std::make_shared<Acts::Transform3D>(
                Acts::Translation3D(pModuleCenter) * pModuleRotation);
            // apply the stereo
            if (m_cfg.posnegModuleBacksideStereo.size()) {
              // twist by the stereo angle
              double stereoBackSide
                  = m_cfg.posnegModuleBacksideStereo.at(ipnl).at(ipnR);
              (*mutableNModuleTransform.get()) *= Acts::AngleAxis3D(
                  -stereoBackSide, Acts::Vector3D::UnitZ());
              (*mutablePModuleTransform.get()) *= Acts::AngleAxis3D(
                  -stereoBackSide, Acts::Vector3D::UnitZ());
            }
            // Finalize the transform
            nModuleTransform = std::const_pointer_cast<const Acts::Transform3D>(
                mutableNModuleTransform);
            pModuleTransform = std::const_pointer_cast<const Acts::Transform3D>(
                mutablePModuleTransform);
            // everything is set for the next module
            FW::Generic::GenericDetectorElement* bsnmodule
                = new FW::Generic::GenericDetectorElement(nModuleIdentifier,
                                                          nModuleTransform,
                                                          moduleBounds,
                                                          moduleThickness,
                                                          moduleMaterialPtr);
            FW::Generic::GenericDetectorElement* bspmodule
                = new FW::Generic::GenericDetectorElement(pModuleIdentifier,
                                                          pModuleTransform,
                                                          moduleBounds,
                                                          moduleThickness,
                                                          moduleMaterialPtr);
            // register the backside of the binmembers
            std::vector<const Acts::DetectorElementBase*> bspbinmember
                = {pmodule};
            std::vector<const Acts::DetectorElementBase*> pbinmember
                = {bspmodule};
            std::vector<const Acts::DetectorElementBase*> bsnbinmember
                = {nmodule};
            std::vector<const Acts::DetectorElementBase*> nbinmember
                = {bsnmodule};
            bsnmodule->registerBinmembers(bsnbinmember);
            nmodule->registerBinmembers(nbinmember);
            bspmodule->registerBinmembers(bspbinmember);
            pmodule->registerBinmembers(pbinmember);
            // memory management - we need a detector store to hold them
            // somewhere @todo add detector store facility
            m_posnegModule.push_back(bsnmodule);
            m_posnegModule.push_back(bspmodule);
          }
          // create the surface
          nsVector.push_back(nmodule->surface().getSharedPtr());
          psVector.push_back(pmodule->surface().getSharedPtr());
        }
        // counter of rings
        ++ipnR;
      }
      // the binning
      size_t layerBinsR = m_cfg.posnegModulePhiBins.at(ipnl).size();
      // never multiply 1 single r-bin, does not make sense
      if (layerBinsR > 1) {
        // multiply with the given bin mulitplier
        layerBinsR *= m_cfg.posnegLayerBinMultipliers.first;
      }
      size_t layerBinsPhi = 0;
      // take the maximum phi bins in that layer
      for (unsigned int phiBins : m_cfg.posnegModulePhiBins.at(ipnl)) {
        layerBinsPhi = phiBins > layerBinsPhi ? phiBins : layerBinsPhi;
        layerBinsPhi *= m_cfg.posnegLayerBinMultipliers.second;
      }
      // create the layers with the surface arrays
      Acts::ProtoLayer pln(nsVector);
      pln.envR = {layerEnvelopeR, layerEnvelopeR};
      pln.envZ = {m_cfg.approachSurfaceEnvelope, m_cfg.approachSurfaceEnvelope};
      Acts::MutableLayerPtr nLayer = m_cfg.layerCreator->discLayer(
          nsVector, layerBinsR, layerBinsPhi, pln);
      Acts::ProtoLayer plp(psVector);
      plp.envR = {layerEnvelopeR, layerEnvelopeR};
      plp.envZ = {m_cfg.approachSurfaceEnvelope, m_cfg.approachSurfaceEnvelope};
      Acts::MutableLayerPtr pLayer = m_cfg.layerCreator->discLayer(
          psVector, layerBinsR, layerBinsPhi, plp);

      // the layer is built le't see if it needs material
      if (m_cfg.posnegLayerMaterialProperties.size()) {
        std::shared_ptr<const Acts::SurfaceMaterial> layerMaterialPtr(
            new Acts::HomogeneousSurfaceMaterial(
                m_cfg.posnegLayerMaterialProperties[ipnl]));
        // central material
        if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) == 0.) {
          // assign the surface material - the layer surface is the material
          // surface
          nLayer->surfaceRepresentation().setAssociatedMaterial(
              layerMaterialPtr);
          pLayer->surfaceRepresentation().setAssociatedMaterial(
              layerMaterialPtr);
          ACTS_VERBOSE("- and material at central layer surface.");
        } else {
          // approach surface material
          // get the approach descriptor - at this stage we know that the
          // approachDescriptor exists
          auto nApproachSurfaces
              = nLayer->approachDescriptor()->containedSurfaces();
          auto pApproachSurfaces
              = pLayer->approachDescriptor()->containedSurfaces();
          if (m_cfg.posnegLayerMaterialConcentration.at(ipnl) > 0.) {
            auto mutableInnerNSurface
                = const_cast<Acts::Surface*>(nApproachSurfaces.at(0));
            mutableInnerNSurface->setAssociatedMaterial(layerMaterialPtr);
            auto mutableOuterPSurface
                = const_cast<Acts::Surface*>(pApproachSurfaces.at(1));
            mutableOuterPSurface->setAssociatedMaterial(layerMaterialPtr);
            ACTS_VERBOSE("- and material at outer approach surfaces.");
          } else {
            auto mutableOuterNSurface
                = const_cast<Acts::Surface*>(nApproachSurfaces.at(1));
            mutableOuterNSurface->setAssociatedMaterial(layerMaterialPtr);
            auto mutableInnerPSurface
                = const_cast<Acts::Surface*>(pApproachSurfaces.at(0));
            mutableInnerPSurface->setAssociatedMaterial(layerMaterialPtr);
            ACTS_VERBOSE("- and material at inner approach surfaces.");
          }
        }
      }
      // push it into the layer vector
      m_nLayers.push_back(nLayer);
      m_pLayers.push_back(pLayer);
    }
  }
}
