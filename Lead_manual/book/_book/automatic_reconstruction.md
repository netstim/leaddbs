## Automatic Selection of Entry Points

By default, _Lead-DBS_ performs an **automatic** search of the artifacts caused by the electrode leads with the following parameters:
- _Entry point:_ STN, GPi, or ViM
- _Axis:_ Use average of coronal and transversal, smoothed
- _Mask window size:_ auto

In most cases, an automatic reconstruction will render an adequate result. If this is not the case, then the user can improve the process by choosing different options within the parameters `Axis`or `Mask window size`.

**Remember:**
Within the parameter `Entry point for electrodes`, only the first two options will run an automatic processing.
