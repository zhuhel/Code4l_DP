################################################################################
# Package: MyAnalysis
################################################################################

SET(Boost_NO_BOOST_CMAKE ON CACHE BOOL "TurnOff BoostFinder" FORCE)

# Declare the package name:
atlas_subdir( MyAnalysis )

# Define a precompiler variable for the release. Rel 20.1 == 2001; 20.7 = 2007; 21.0 - 2100
add_definitions( -D__RELEASE__=2100 )

# Declare the package's dependencies:
atlas_depends_on_subdirs( PUBLIC
	PhysicsAnalysis/D3PDTools/EventLoopGrid
	PhysicsAnalysis/D3PDTools/EventLoop
	PhysicsAnalysis/D3PDTools/EventLoopAlgs
	Event/xAOD/xAODEventInfo
	Event/xAOD/xAODJet
	Control/xAODRootAccess
	Event/xAOD/xAODMuon
	Event/xAOD/xAODTracking
	Control/AthLinks
	Event/xAOD/xAODCaloEvent
	Control/AthContainers
	Event/xAOD/xAODEgamma
	Event/xAOD/xAODMissingET
	PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonMomentumCorrections
	Event/xAOD/xAODTruth
	PhysicsAnalysis/MuonID/MuonSelectorTools
	PhysicsAnalysis/AnalysisCommon/FsrUtils
	Trigger/TrigAnalysis/TrigDecisionTool
	Event/xAOD/xAODTrigMuon
	Trigger/TrigConfiguration/TrigConfxAOD
	PhysicsAnalysis/AnalysisCommon/IsolationSelection
	DataQuality/GoodRunsLists
	PhysicsAnalysis/AnalysisCommon/PileupReweighting
	Event/xAOD/xAODBTagging
	AssociationUtils/AssociationUtils
	Event/xAOD/xAODCutFlow
	Event/xAOD/xAODCore
	Tools/PathResolver
        #IsolationTool/IsolationTool
        #RecoToolInterfaces/RecoToolInterfaces
        Reconstruction/MET/METUtilities
        PhysicsAnalysis/HiggsPhys/Run2/HZZ/Tools/ZMassConstraint
        PhysicsAnalysis/TopPhys/QuickAna
        #QuickAna/QuickAna
)


# Find the needed external(s):
find_package( ROOT COMPONENTS Core RIO Hist Tree )

# build a CINT dictionary for the library
atlas_add_root_dictionary ( MyAnalysisLib MyAnalysisCintDict
	ROOT_HEADERS MyAnalysis/MyxAODAnalysis.h Root/LinkDef.h
                            EXTERNAL_PACKAGES ROOT
)

# build a shared library
atlas_add_library( MyAnalysisLib
	MyAnalysis/*.h Root/*  ${MyAnalysisCintDict}
                   PUBLIC_HEADERS MyAnalysis
                   INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   LINK_LIBRARIES ${ROOT_LIBRARIES}
                   EventLoop
                   EventLoopAlgs
                   EventLoopGrid
                   GoodRunsListsLib
                   IsolationSelectionLib
                   MuonEfficiencyCorrectionsLib
                   MuonMomentumCorrectionsLib
                   MuonSelectorToolsLib
                   PathResolver
                   PileupReweightingLib
                   TrigConfxAODLib
                   TrigDecisionToolLib
                   xAODBTagging
                   xAODEgamma
                   xAODEventInfo
                   xAODJet
                   xAODMissingET
                   xAODMuon
                   xAODRootAccess
                   xAODTracking
                   xAODTrigEgamma
                   xAODTruth
		   xAODCutFlow
		   FsrUtilsLib
		   QuickAnaLib
		   ElectronPhotonSelectorToolsLib
		   ElectronPhotonFourMomentumCorrectionLib
		   TriggerMatchingToolLib
                   METUtilitiesLib
                   #RecoToolInterfacesLib
                   AssociationUtilsLib
                   ZMassConstraintLib
)

if( XAOD_STANDALONE )
   foreach( exec runFourLep
            grid_submission
            grid_data16 grid_data15
            grid_data17_periodB grid_data17_periodC grid_data17_periodD 
            grid_data17_periodE grid_data17_periodF grid_data17_periodH 
            grid_data17_periodI grid_data17_periodK )
      atlas_add_executable( ${exec}
         util/${exec}.cxx
         INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}
         LINK_LIBRARIES ${ROOT_LIBRARIES} ${Boost_LIBRARIES}
         MyAnalysisLib)
   endforeach()
endif()

# Install files from the package:
#atlas_install_scripts( scripts/*.py )
atlas_install_data( share/* )
