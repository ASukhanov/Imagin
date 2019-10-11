# Imagin
Interactive Viewer/Analyzer for streamed images or files based on python, pyqtgraph and scipy.

Features:
+ Input source: 
    * File system, 
    * HTTP link to image, 
    * RHIC ADO parameter (requires cns.py), 
    * EPICS PV (requires pyepics),
    * Local USB camera (requires openCV).
+ Default image format is PNG, other formats supported as well.
+ 16+ bit/channel images supported.
+ Image orientation and rotation (program options: -o and -R).
+ Interactive zooming, panning, rotation.
+ Contrast/Coloration control.
+ Region Of Interest (ROI) for image analysis.
+ ROI partitioning. 
+ Isocurve. The isocurve level defines the threshold for object finding.
+ Fast and robust characterization of multiple objects using fitted ellipsoids.
+ Gaussian (1D and 2D) fit for improved precision.
+ Interactive python console with access to image data, graphics objects and shell.
+ Interactive calibration of pixels to millimeters.
+ Reference Images: save/retrieve image to/from a reference slots.
+ Background subtraction using a reference image.
+ Easy extentable using with user-suplied add-ons.
+ Fast browsing/cleanup of image directories.
## Try:
    python imagin.py
![](imagin_screenshot.png)

