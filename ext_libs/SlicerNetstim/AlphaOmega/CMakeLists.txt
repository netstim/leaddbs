# BUILD MODULE https://www.slicer.org/wiki/Documentation/Nightly/Developers/Build_Module
# cmake_minimum_required(VERSION 2.8.9)
# find_package(Slicer COMPONENTS ConfigurePrerequisites) 
# project(EmptyExtensionTemplate)

#-----------------------------------------------------------------------------
set(MODULE_NAME AlphaOmega)
set(MODULE_TITLE ${MODULE_NAME})

string(TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER)

# BUILD MODULE https://www.slicer.org/wiki/Documentation/Nightly/Developers/Build_Module
# find_package(Slicer REQUIRED) 
# include(${Slicer_USE_FILE}) 

#-----------------------------------------------------------------------------
# Alpha Omega SDK

SET(ALPHAOMEGA_SDK_PATH_HINTS
  "C:/Program Files (x86)/AlphaOmega/Neuro Omega System SDK"
  )
find_path(ALPHAOMEGA_SDK_DIR ReadMe.txt
  PATHS ${ALPHAOMEGA_SDK_PATH_HINTS}
  DOC "AlphaOmega SDK directory"
  )

#-----------------------------------------------------------------------------
add_subdirectory(MRML)
add_subdirectory(Logic)

#-----------------------------------------------------------------------------
set(MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT")

# Current_{source,binary} and Slicer_{Libs,Base} already included
set(MODULE_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/Logic
  ${CMAKE_CURRENT_BINARY_DIR}/Logic
  ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/Include
  )

  
#-----------------------------------------------------------------------------

set(MODULE_SRCS
  qSlicer${MODULE_NAME}Module.cxx
  qSlicer${MODULE_NAME}Module.h
  qSlicer${MODULE_NAME}ModuleWidget.cxx
  qSlicer${MODULE_NAME}ModuleWidget.h
  qAlphaOmegaStatusThread.cxx
  qAlphaOmegaStatusThread.h
  )

set(MODULE_MOC_SRCS
  qSlicer${MODULE_NAME}Module.h
  qSlicer${MODULE_NAME}ModuleWidget.h
  qAlphaOmegaStatusThread.h
  )

set(MODULE_UI_SRCS
  Resources/UI/qSlicer${MODULE_NAME}ModuleWidget.ui
  )

set(MODULE_TARGET_LIBRARIES
  vtkSlicer${MODULE_NAME}ModuleLogic
  ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64/NeuroOmega_x64.lib
  )

set(MODULE_RESOURCES
  Resources/qSlicer${MODULE_NAME}Module.qrc
  )

#-----------------------------------------------------------------------------
slicerMacroBuildLoadableModule(
  NAME ${MODULE_NAME}
  TITLE ${MODULE_TITLE}
  EXPORT_DIRECTIVE ${MODULE_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  SRCS ${MODULE_SRCS}
  MOC_SRCS ${MODULE_MOC_SRCS}
  UI_SRCS ${MODULE_UI_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RESOURCES ${MODULE_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
add_custom_command(
  TARGET qSlicer${MODULE_NAME}Module POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64/NeuroOmega_x64.lib
    ${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_BIN_DIR}/Release/NeuroOmega_x64.lib)

add_custom_command(
  TARGET qSlicer${MODULE_NAME}Module POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64/NeuroOmega_x64.dll
    ${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_BIN_DIR}/Release/NeuroOmega_x64.dll)