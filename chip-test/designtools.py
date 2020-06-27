"""
    Hi👌
"""
import numpy as np

def gen_wafer(radius=1000, delta=np.pi/10, flip_wafer=True):
    """Generete wafer semi-circle"""
    angles = np.linspace(delta, 2*np.pi-delta, 81)
    if flip_wafer:
        angles += np.pi
    xs = radius * np.cos(angles)
    ys = radius * np.sin(angles)
    return list(zip(xs,ys))

def dice_mark(x, y, angle=0, flip_x=False, flip_y=False):
    """Generate single dice mark (staircase thing)"""
    pts = np.array([
        [0, 0], [300, 0], [300, 50], [250, 50], [250, 60], [200, 60], [200, 70],
        [150, 70], [150, 80], [100, 80], [100, 90], [50, 90], [50, 100], [0, 100]
    ])
    pts[:, 0] -= 300
    pts[:, 1] -= 100
    # Flip
    pts[:, 0] = -pts[:, 0] if flip_x else pts[:, 0]
    pts[:, 1] = -pts[:, 1] if flip_y else pts[:, 1]
    # Rotate
    pts = pts@rot_mat(angle)
    # Move
    pts[:, 0] += x
    pts[:, 1] += y
    return pts

def tri(x, y, dx, dy):
    """Create a triangle"""
    return [[x-dx/2, y+dy], [x+dx/2, y+dy], [x, y]]

def rot_mat(angle):
    """Return the 2D rotation matrix"""
    return np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle), np.cos(angle)]])
