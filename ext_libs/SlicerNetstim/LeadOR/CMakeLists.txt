# -----------------------------------------------------------------------------
set(MODULE_NAME LeadOR)

# -----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  LeadORLib/Widgets/tables.py
  LeadORLib/Widgets/__init__.py
  LeadORLib/util.py
  LeadORLib/__init__.py
)

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/Icons/LeadOR1.png
  Resources/Icons/VolumeResliceDriver.png
  Resources/Icons/VcrRecord16.png
  Resources/UI/${MODULE_NAME}.ui
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
