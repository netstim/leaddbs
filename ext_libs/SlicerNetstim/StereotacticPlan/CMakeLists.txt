# -----------------------------------------------------------------------------
set(MODULE_NAME StereotacticPlan)

# -----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  StereotacticPlanLib/ImportFrom/Import_From_Brainlab.py
  StereotacticPlanLib/ImportFrom/Import_From_ROSA.py
  StereotacticPlanLib/ImportFrom/importerBase.py
  StereotacticPlanLib/ImportFrom/__init__.py
  StereotacticPlanLib/Widgets/CustomWidgets.py
  StereotacticPlanLib/Widgets/__init__.py
)

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/UI/${MODULE_NAME}.ui
  Resources/Icons/VolumeResliceDriver.png
)

# -----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
)

# -----------------------------------------------------------------------------
if(BUILD_TESTING)
  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
