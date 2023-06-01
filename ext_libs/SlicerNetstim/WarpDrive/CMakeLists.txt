#-----------------------------------------------------------------------------
set(MODULE_NAME WarpDrive)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  WarpDriveLib/Effects/CircleEffect.py
  WarpDriveLib/Effects/DrawEffect.py
  WarpDriveLib/Effects/Effect.py
  WarpDriveLib/Effects/PointerEffect.py
  WarpDriveLib/Effects/PointToPointEffect.py
  WarpDriveLib/Effects/__init__.py
  WarpDriveLib/Helpers/GridNodeHelper.py
  WarpDriveLib/Helpers/LeadDBSCall.py
  WarpDriveLib/Helpers/__init__.py
  WarpDriveLib/Tools/DrawTool.py
  WarpDriveLib/Tools/NoneTool.py
  WarpDriveLib/Tools/PointToPointTool.py
  WarpDriveLib/Tools/SmudgeTool.py
  WarpDriveLib/Tools/__init__.py
  WarpDriveLib/Widgets/Buttons.py
  WarpDriveLib/Widgets/Toolbar.py
  WarpDriveLib/Widgets/ToolWidget.py
  WarpDriveLib/Widgets/Tables.py
  WarpDriveLib/Widgets/__init__.py
  WarpDriveLib/__init__.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/Icons/WarpDrive1.png
  Resources/Icons/Add.png
  Resources/Icons/Delete.png
  Resources/Icons/Draw.png
  Resources/Icons/None.png
  Resources/Icons/PointToPoint.png
  Resources/Icons/Rename.png
  Resources/Icons/Smudge.png
  Resources/Icons/SlicerVisible.png
  Resources/Icons/SlicerUndo.png
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
