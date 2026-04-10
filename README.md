# Imagin. (Beta release).

Imagin is an interactive viewer and analyzer for streamed images and image files, built with Python, pyqtgraph, and SciPy.

## Features

- Input sources:
  - file system
  - HTTP image URLs
  - EPICS Channel Access
  - EPICS PVAccess
  - USB cameras
- PNG is the default image format; other formats are also supported.
- Supports images with 16-bit or higher channel depth.
- Image orientation control and arbitrary rotation.
- Interactive zooming and panning.
- Contrast and color-map control.
- Region of interest (ROI) selection for image analysis.
- ROI partitioning.
- ROI projection plots.
- Isocurves, with the isocurve level used as the threshold for object detection.
- Fast, robust characterization of multiple objects using fitted ellipsoids.
- Optional Gaussian fitting in 1D or 2D for improved precision.
- De-speckling.
- Automatic background elimination using prominence filtering (de-base).
- Background subtraction using reference images.
- Interactive calibration from pixels to millimeters.
- Reference image slots for saving and restoring images.
- Extensible through user-supplied add-ons.
- Fast browsing and cleanup of image directories.
- Interactive Python console with access to image data, graphics objects, and the shell.

![Imagin screenshot](docs/imagin_screenshot.png)

[Presentation](docs/Slides_from_ICALEPCS-2019.pdf)

## Examples

```bash
python -m imagin -b file sample_images/*.jpg -t 100 -m 40
python -m imagin -b file ~/Pictures/*.png
python -m imagin -b pva image0:image
```
To start the image simulator for `image0:image`:
```bash
pip install epicsdev
python -m epicsdev.imagegen
```

The simulated EPICS camera can also be run from Docker:

https://hub.docker.com/r/klauer/simioc-docker

```bash
python3 setup_SimDetector.py
python3 -m imagin -b epics sim:det -m 16 -t 80
```

