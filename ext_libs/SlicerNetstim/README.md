This branch contains only the scripted modules of the extension. 

It is used to copy the modules to the Lead-DBS ext_libs folder:

`git subtree pull --prefix ext_libs/SlicerNetstim --squash https://github.com/netstim/SlicerNetstim.git scripted_modules_only`

They are loaded when running the [custom Slicer app for Lead-DBS](https://github.com/netstim/SlicerForLeadDBS) from Lead-DBS.

This allows to update the scripted modules easier than compiling the custom Slicer every time they need to be updated.