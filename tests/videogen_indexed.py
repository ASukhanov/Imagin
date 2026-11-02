import cv2 as cv
import numpy as np

scale = 10
WIDTH, HEIGHT = int(100 * scale), int(100 * scale)
FPS = 10
NUM_FRAMES = 100

rng = np.random.default_rng(seed=42)

# Lossless 16-bit grayscale.
ffv1 = cv.VideoWriter(
    "random_ffv1.mkv",
    cv.CAP_FFMPEG,
    cv.VideoWriter_fourcc(*"FFV1"),
    FPS,
    (WIDTH, HEIGHT),
    [
        cv.VIDEOWRITER_PROP_DEPTH, cv.CV_16U,
        cv.VIDEOWRITER_PROP_IS_COLOR, 0,
    ],
)

# Lossy 8-bit grayscale.
mp4 = cv.VideoWriter(
    "random_mp4.mp4",
    cv.CAP_FFMPEG,
    cv.VideoWriter_fourcc(*"mp4v"),
    FPS,
    (WIDTH, HEIGHT),
    False,  # isColor
)

FONT = cv.FONT_HERSHEY_SIMPLEX
THICKNESS = int(2*scale)
TEXT_HEIGHT = int(20*scale)
PADDING = int(5*scale)

font_scale = cv.getFontScaleFromHeight(
    FONT, TEXT_HEIGHT, THICKNESS
)

try:
    if not ffv1.isOpened():
        raise RuntimeError(
            "Cannot open FFV1 writer; requires OpenCV with "
            "FFmpeg and 16-bit video writing support."
        )
    if not mp4.isOpened():
        raise RuntimeError("Cannot open MPEG-4 writer.")

    for i in range(NUM_FRAMES):
        frame16 = rng.integers(
            0, 65536,
            size=(HEIGHT, WIDTH),
            dtype=np.uint16,
        )

        text = str(i + 1)
        (text_width, text_height), baseline = cv.getTextSize(
            text, FONT, font_scale, THICKNESS
        )

        # # Black background behind the number.
        # frame16[
        #     :text_height + baseline + 2 * PADDING,
        #     :text_width + 2 * PADDING,
        # ] = 0

        # Render bold, 20-pixel-high text into an 8-bit mask.
        mask = np.zeros((HEIGHT, WIDTH), dtype=np.uint8)
        cv.putText(
            mask,
            text,
            (PADDING, PADDING + text_height),
            FONT,
            font_scale,
            color=255,
            thickness=THICKNESS,
            lineType=cv.LINE_8,
        )

        # Transfer white text to the 16-bit frame.
        frame16[mask != 0] = 65535

        ffv1.write(frame16)

        frame8 = (frame16 >> 8).astype(np.uint8)
        mp4.write(frame8)

finally:
    ffv1.release()
    mp4.release()

from pathlib import Path

for name in ("random_ffv1.mkv", "random_mp4.mp4"):
    print(f"{name}: {Path(name).stat().st_size / 1024:.1f} KiB")

