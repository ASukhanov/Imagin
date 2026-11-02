"""Frame-by-frame video viewer for 16-bit and 8-bit videos.
 It uses OpenCV to read and display video files, allowing the user to navigate
 through frames using keyboard input. Once it reaches the end of the video, it
 will loop back to the beginning and continue displaying frames in reverse order.
"""
__version__ = 'v0.0.2 2026-10-26'# 
import os
os.environ["QT_QPA_PLATFORM"] = "xcb"# to avoid "Could not connect to any X display" error when running in headless mode
import argparse

import cv2 as cv
import numpy as np

def display_video(filename):
    """Display a video frame by frame.

    Args:
        filename (str): Path to the video file.
        is_16bit (bool): Whether the video is 16-bit grayscale.

    Returns:
        bool: True if the video was displayed completely, False if interrupted.
    """
    cap = cv.VideoCapture(filename, cv.CAP_FFMPEG)
    window = "Frame viewer"

    try:
        if not cap.isOpened():
            raise RuntimeError(f"Cannot open {filename}")

        # if is_16bit:
        #     # Disable automatic conversion to 8-bit BGR.
        #     cap.set(cv.CAP_PROP_CONVERT_RGB, 0)

        cv.namedWindow(window, cv.WINDOW_NORMAL)
        cv.resizeWindow(window, 500, 500)

        frame_count = int(cap.get(cv.CAP_PROP_FRAME_COUNT))
        increment = 1
        index = -1
        while True:
            index += increment
            if index < 0:
                increment = 1
                index = 0
            elif index >= frame_count:
                increment = -1
                index = frame_count - 1
            cap.set(cv.CAP_PROP_POS_FRAMES, index)
            ok, frame = cap.read()
            if not ok:
                break

            if frame.dtype != np.uint16 or frame.ndim != 2:
                 display = frame
            else:
                # Original 16-bit data remains available in `frame`.
                display = (frame >> 8).astype(np.uint8)

            cv.setWindowTitle(
                window, f"{filename} — frame {index}/{frame_count-1}"
            )
            cv.imshow(window, display)
            key = cv.waitKey(0) & 0xFF
            if key in (ord("q"), ord("Q"), 27):
                return False

        return True

    finally:
        cap.release()
        cv.destroyAllWindows()


#if display_video("random_ffv1.mkv", is_16bit=True):
##    display_video("random_mp4.mp4")

parser = argparse.ArgumentParser(description=__doc__
,formatter_class=argparse.ArgumentDefaultsHelpFormatter
,epilog=(f'{__version__},'))

parser.add_argument("filename", type=str, help=
                    "Path to the video file.",)
args = parser.parse_args()

display_video(args.filename)
