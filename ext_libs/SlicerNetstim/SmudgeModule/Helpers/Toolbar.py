import qt
import slicer


class reducedToolbar(object):

  def __init__(self):

    self.toolbar = qt.QToolBar()

    #self.setFixedHeight(21)

    fourUpPixmap = qt.QPixmap('/Users/simon/Desktop/Icons/LayoutFourUpView.png')
    fourUpIcon = qt.QIcon(fourUpPixmap)
    self.layoutButton = qt.QToolButton(self.toolbar)
    self.layoutButton.setIcon(fourUpIcon)
    self.layoutButton.setFixedSize(fourUpPixmap.rect().size())
    self.layoutButton.setPopupMode(2)

    oneUpPixmap = qt.QPixmap('/Users/simon/Desktop/Icons/LayoutOneUpView.png')
    oneUpIcon = qt.QIcon(oneUpPixmap)

    menu = qt.QMenu()
    allSlicesAction = menu.addAction(fourUpIcon, 'All Slices')
    axialAction     = menu.addAction(oneUpIcon,  'Axial')
    sagittalAction  = menu.addAction(oneUpIcon,  'Sagittal')
    coronalAction   = menu.addAction(oneUpIcon,  'Coronal')

    self.layoutButton.setMenu(menu)

    allSlicesAction.connect('triggered(bool)', self.allSlicesTriggered)


  def getToolbar(self):
    return self.toolbar

  def allSlicesTriggered(self, t):
    layoutManager = slicer.app.layoutManager()
    layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)


    #
    # Custom Toolbar Area
    #

    tb = qt.QToolBar()
    fourUpPixmap = qt.QPixmap('/Users/simon/Desktop/Icons/LayoutFourUpView.png')
    fourUpIcon = qt.QIcon(fourUpPixmap)
    layoutButton = qt.QToolButton(tb)
    layoutButton.setIcon(fourUpIcon)
    layoutButton.setFixedWidth(23)
    layoutButton.setPopupMode(2)

    oneUpPixmap = qt.QPixmap('/Users/simon/Desktop/Icons/LayoutOneUpView.png')
    oneUpIcon = qt.QIcon(oneUpPixmap)

    menu = qt.QMenu()
    allSlicesAction = menu.addAction(fourUpIcon, 'All Slices')
    axialAction     = menu.addAction(oneUpIcon,  'Axial')
    sagittalAction  = menu.addAction(oneUpIcon,  'Sagittal')
    coronalAction   = menu.addAction(oneUpIcon,  'Coronal')
    layoutButton.setMenu(menu)

    slicer.util.mainWindow().addToolBar(tb)