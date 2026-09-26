import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from PIL import Image
import os

# ------------------------------------------------------------
# Abflachung der Sonne – Python-Portierung von AbflachungSonne.m
# ------------------------------------------------------------

# Einmalige Parametrierung des Bildpfads
IMG_PATH = "../../hm/Data/"   # <--- hier liegen SunsetHawaii.jpg usw.

# Farben ähnlich MATLAB GetColorLines
def get_color_lines(n=10):
    cmap = plt.get_cmap("tab10")
    return np.array([cmap(i) for i in range(n)])

Colors = get_color_lines(12)

# Sun elevation values
h = np.linspace(1.25, -1.75, 6)  # degrees

# Sun outline
phi = np.linspace(0, 2*np.pi, 200)
Dsun = 32.3/60.0  # Sun radius in degrees
cphi = 0.5 * Dsun * np.cos(phi)
sphi = 0.5 * Dsun * np.sin(phi)

# Load picture(s)
pic_flag = 0

if pic_flag:
    rgb1 = np.array(Image.open(os.path.join(IMG_PATH, "SunsetHawaii.jpg"))).astype(float)
    rgb2 = np.array(Image.open(os.path.join(IMG_PATH, "SunsetHawaii01.jpg"))).astype(float)
    rgb = (rgb1 + rgb2) / 2.0
else:
    rgb = np.array(Image.open(os.path.join(IMG_PATH, "SunsetHawaii.jpg"))).astype(float)

pixres = 0.00101  # degrees per pixel
xra = pixres * rgb.shape[0] / 2
yra = pixres * rgb.shape[1] / 2

# Flip image vertically (MATLAB flip(rgb,1))
rgb_flip = np.flip(rgb, axis=0)

fig, ax = plt.subplots(figsize=(8, 8))

ax.imshow(
    rgb_flip,
    extent=[-xra, xra, -yra + 0.605, yra + 0.605],
    alpha=0.6,
    origin="upper"
)

ax.grid(True)
sh = 0.6  # shift for plot

# Plot sun outlines for different elevations
for k in range(len(h)):
    cys = h[k] + sphi
    R1 = 1.02 / np.tan(np.deg2rad(cys + 10.3 / (cys + 5.11))) / 60.0
    cy = cys + R1

    if k < len(h) - 1:
        ax.plot(cphi - sh, cy, color=Colors[k], linewidth=3, linestyle=":")
        poly = Polygon(np.column_stack((cphi - sh, cy)),
                       closed=False,
                       facecolor=Colors[10],
                       alpha=0.75)
        ax.add_patch(poly)

    ax.plot(cphi + sh, cys, color=Colors[k], linewidth=2)

# Location example
h0 = -0.3492
cys0 = h0 + sphi
R1 = 1.02 / np.tan(np.deg2rad(cys0 + 10.3 / (cys0 + 5.11))) / 60.0
cy0 = cys0 + R1
ax.plot(cphi, cy0, color="k", linewidth=2, linestyle=":")

# Horizon line
hor = -0.0907
ax.plot([-1, 1], [hor, hor], color=Colors[3], linewidth=2)

xh = [-1, 1, 1, -1]
yh = [-1, -1, hor, hor]
horizon_patch = Polygon(np.column_stack((xh, yh)),
                        closed=True,
                        facecolor=Colors[8],
                        alpha=0.65)
ax.add_patch(horizon_patch)

ax.set_aspect("equal")
ax.set_ylabel("Höhe in °", fontsize=16)
ax.set_xlabel("Azimut (relativ) in °", fontsize=16)
ax.set_xlim([-0.88, 0.88])
ax.set_ylim([-0.27, 1.49])
ax.tick_params(labelsize=14)

plt.show()

