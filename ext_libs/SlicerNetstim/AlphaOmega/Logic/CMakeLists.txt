project(vtkSlicer${MODULE_NAME}ModuleLogic)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  ${MRML_LIBRARIES}
  ${vtkSlicer${MODULE_NAME}ModuleMRML_SOURCE_DIR}
  ${vtkSlicer${MODULE_NAME}ModuleMRML_BINARY_DIR}
  ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/Include
  )

set(${KIT}_SRCS
  vtkSlicer${MODULE_NAME}Logic.cxx
  vtkSlicer${MODULE_NAME}Logic.h
  )

set(${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${MRML_LIBRARIES}
  vtkSlicer${MODULE_NAME}ModuleMRML
  ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64/NeuroOmega_x64.lib
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  )



# if (ALPHAOMEGA_SDK_DIR)
#   include_directories(${ALPHAOMEGA_SDK_DIR}/CPP_SDK/Include)
#   link_directories(${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64) # TODO: x32
#   target_link_libraries(vtkSlicer${MODULE_NAME}ModuleLogic ${ALPHAOMEGA_SDK_DIR}/CPP_SDK/win64/NeuroOmega_x64.lib)
# endif()

