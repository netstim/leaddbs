import slicer
import os

extensionName = 'SlicerNetstim'
em = slicer.app.extensionsManagerModel()

print('os:%s' % em.slicerOs)
print('rev:%s' % em.slicerRevision)
print('extensionsInstallPath:%s' % em.extensionsInstallPath())
print('slicerMajorVersion:%s' % slicer.app.majorVersion)

# exit with error if takes longer than a minute
qt.QTimer.singleShot(60*1000, lambda: slicer.util.exit(1))

if not em.isExtensionInstalled(extensionName):
    extensionMetaData = em.retrieveExtensionMetadataByName(extensionName)
    if len(extensionMetaData) == 0:
        slicer.util.exit(1)
    url = f"{em.serverUrl().toString()}/api/v1/item/{extensionMetaData['_id']}/download"
    extensionPackageFilename = os.path.join(slicer.app.temporaryPath, extensionMetaData['_id'])
    slicer.util.downloadFile(url, extensionPackageFilename)
    em.interactive = False  # Disable popups (automatically install dependencies)
    em.installExtension(extensionPackageFilename)

slicer.util.exit()