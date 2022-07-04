#-----------------------------------------------------------------------------
set(MODULE_NAME WarpDrive)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  WarpDriveLibLocal/Effects/CircleEffect.py
  WarpDriveLibLocal/Effects/DrawEffect.py
  WarpDriveLibLocal/Effects/Effect.py
  WarpDriveLibLocal/Effects/PointerEffect.py
  WarpDriveLibLocal/Effects/PointToPointEffect.py
  WarpDriveLibLocal/Effects/__init__.py
  WarpDriveLibLocal/Helpers/GridNodeHelper.py
  WarpDriveLibLocal/Helpers/LeadDBSCall.py
  WarpDriveLibLocal/Helpers/__init__.py
  WarpDriveLibLocal/Tools/DrawTool.py
  WarpDriveLibLocal/Tools/NoneTool.py
  WarpDriveLibLocal/Tools/PointToPointTool.py
  WarpDriveLibLocal/Tools/SmudgeTool.py
  WarpDriveLibLocal/Tools/__init__.py
  WarpDriveLibLocal/Widgets/Toolbar.py
  WarpDriveLibLocal/Widgets/ToolWidget.py
  WarpDriveLibLocal/Widgets/Tables.py
  WarpDriveLibLocal/Widgets/__init__.py
  WarpDriveLibLocal/__init__.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/Icons/Add.png
  Resources/Icons/Delete.png
  Resources/Icons/Draw.png
  Resources/Icons/None.png
  Resources/Icons/PointToPoint.png
  Resources/Icons/Rename.png
  Resources/Icons/Smudge.png
  Resources/Icons/SlicerVisible.png
  Resources/UI/${MODULE_NAME}.ui
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
