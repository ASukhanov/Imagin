# Imagin
Interactive Viewer/Analyzer for streamed images or files based on python, pyqtgraph and scipy.

Features:
+ Input source: 
    * File system, 
    * HTTP link to image,
    * EPICS PV,
    * Local USB cameras (requires openCV).
+ Default image format is PNG, other formats supported as well.
+ 16+ bit/channel images supported.
+ Image orientation and arbitrary rotation.
+ Interactive zooming, panning, rotation.
+ Contrast/Coloration control.
+ Region Of Interest (ROI) for image analysis.
+ ROI partitioning. 
+ ROI projection plots.
+ Isocurves. The isocurve level defines the threshold for object finding.
+ Fast and robust characterization of multiple objects using fitted ellipsoids.
+ Gaussian (1D and 2D) fit for improved precision.
+ Effective de-speckling.
+ Automatic elimination of background using prominence filter (de-base).
+ Background subtraction using a reference image.
+ Interactive calibration of pixels to millimeters.
+ Reference Images: save/retrieve image to/from a reference slots.
+ Easy extentable using with user-suplied add-ons.
+ Fast browsing/cleanup of image directories.
+ Interactive python console with access to image data, graphics objects and shell.

![](docs/imagin_screenshot.png)

[Presentation](docs/Slides_from_ICALEPCS-2019.pdf)
## Examples

    python3 -m imagin -bfile docs/GalaxyClusterAbell1689_sn_l.jpg -t100 -m40

Note: the simulated camera under EPICS can be run from a docker:
https://hub.docker.com/r/klauer/simioc-docker

Set up the simulated camera:

    python3 setup_SimDetector.py

Start application:

    python3 -m imagin -b epics sim:det -m16 -t80

