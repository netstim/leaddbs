import qt

class ImporterDialogBase():
    def __init__(self):
        self.importerName = None
        self.fileSelectTitle = None
        self.fileSelectExt = None

        self.selectedFile = None
        self.computeReferenceToFrame = False
        self.importACPCCoordinates = False
        self.importDICOM = False
        self.DICOMDir = False

    def run(self, importInFrameSpace):
        dialog = qt.QDialog()
        dialog.setWindowTitle(self.importerName + ' Import Options')

        fileSelectButton = qt.QPushButton('Click to select')
        fileSelectButton.clicked.connect(lambda: fileSelectButton.setText(qt.QFileDialog.getOpenFileName(qt.QWidget(), self.fileSelectTitle, '', '*.' + self.fileSelectExt)))

        computeReferenceToFrameCheckBox = qt.QCheckBox()
        computeReferenceToFrameCheckBox.setEnabled(False)

        importACPCCoordinatesCheckBox = qt.QCheckBox()
        importACPCCoordinatesCheckBox.connect("toggled(bool)", lambda b: computeReferenceToFrameCheckBox.setEnabled(b))

        DICOMDirButton = qt.QPushButton('Click to select')
        DICOMDirButton.clicked.connect(lambda: DICOMDirButton.setText(qt.QFileDialog.getExistingDirectory(qt.QWidget(), 'Select DICOM directory', '')))
        DICOMDirButton.setEnabled(False)

        importDICOMCheckBox = qt.QCheckBox()
        importDICOMCheckBox.connect("toggled(bool)", lambda b: DICOMDirButton.setEnabled(b))

        buttonBox = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel, qt.Qt.Horizontal, dialog)
        buttonBox.accepted.connect(lambda: dialog.accept())
        buttonBox.rejected.connect(lambda: dialog.reject())

        form = qt.QFormLayout(dialog)
        form.addRow('Select File: ', fileSelectButton)
        if not importInFrameSpace:
            form.addRow('Import ACPC coords: ', importACPCCoordinatesCheckBox)
            form.addRow('Compute reference to frame transform: ', computeReferenceToFrameCheckBox)
            form.addRow('Import reference image: ', importDICOMCheckBox)
            form.addRow('DICOM directory: ', DICOMDirButton)
        form.addRow(buttonBox)

        if dialog.exec() == qt.QDialog.Accepted:
            self.selectedFile = fileSelectButton.text
            self.computeReferenceToFrame = computeReferenceToFrameCheckBox.checked
            self.importACPCCoordinates = importACPCCoordinatesCheckBox.checked
            self.importDICOM = importDICOMCheckBox.checked
            self.DICOMDir = DICOMDirButton.text
