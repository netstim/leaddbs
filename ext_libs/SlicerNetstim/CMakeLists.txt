cmake_minimum_required(VERSION 3.13.4)

project(SlicerNetstim)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/SlicerNetstim")
set(EXTENSION_CATEGORY "Netstim")
set(EXTENSION_CONTRIBUTORS "Simon Oxenford (Charite Berlin)")
set(EXTENSION_DESCRIPTION "Netstim modules collection")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SlicerNetstim.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/SlicerNetstim/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(LeadOR)
#add_subdirectory(AlphaOmega) # don't add for slicer
add_subdirectory(StereotacticPlan)
add_subdirectory(ImportAtlas)
add_subdirectory(WarpDrive)
add_subdirectory(NetstimPreferences)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
