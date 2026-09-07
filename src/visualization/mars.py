"""Mars disk/sphere and axisymmetric empirical BS/MPB, coordinates in Rm."""
import numpy as np
from matplotlib.patches import Circle

# Same conic parameters as py_space_zc.maven.bs_mpb; no runtime dependency.
MODELS = {"BS": ((0.6, 2.081, 1.026),),
          "MPB": ((0.640, 1.080, 0.770), (1.600, 0.528, 1.009))}


def boundary_profile(name, xmin=-5., samples=500):
    """Return separate (x,rho) branches in Rm. Exclude unphysical r<=0."""
    if name not in MODELS or samples < 3 or not np.isfinite(xmin):
        raise ValueError("Expected BS/MPB, finite xmin and samples >= 3")
    branches = []
    theta = np.linspace(0, np.pi, samples)
    for k, (x0, length, eccentricity) in enumerate(MODELS[name]):
        denominator = 1 + eccentricity*np.cos(theta)
        valid = denominator > 0
        r = length / denominator[valid]
        x, rho = x0+r*np.cos(theta[valid]), r*np.sin(theta[valid])
        keep = x >= xmin
        if name == "MPB":
            keep &= x >= 0 if k == 0 else x < 0
        branches.append((x[keep], rho[keep]))
    return branches


def plot_mars(ax, radius=1., center=None, color="#b97956", alpha=1.,
              texture=False, texture_path=None, facecolor=None, edgecolor="#654b3c",
              lw=.6, zorder=5):
    """Draw a solid disk or sphere. Lengths use the axes' units (normally Rm)."""
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("radius must be positive")
    color = facecolor or color
    if texture:
        if texture_path is None or ax.name == "3d":
            raise ValueError("Texture needs a supplied disk image and a 2D axis")
        import matplotlib.pyplot as plt
        cx, cy = (0,0) if center is None else center
        return ax.imshow(plt.imread(texture_path), extent=(cx-radius,cx+radius,cy-radius,cy+radius), alpha=alpha,zorder=zorder)
    if ax.name == "3d":
        center = np.zeros(3) if center is None else np.asarray(center)
        u, v = np.meshgrid(np.linspace(0, 2*np.pi, 60), np.linspace(0, np.pi, 30))
        return ax.plot_surface(center[0]+radius*np.cos(u)*np.sin(v),
            center[1]+radius*np.sin(u)*np.sin(v), center[2]+radius*np.cos(v),
            color=color, alpha=alpha, linewidth=0, shade=True)
    center = (0, 0) if center is None else center
    return ax.add_patch(Circle(center, radius, facecolor=color, edgecolor=edgecolor,
                               lw=lw, alpha=alpha, zorder=zorder))


def bs_mpb(ax, plane="XZ", draw_bs=True, draw_mpb=True, xmin=-5.,
           x_slice=None, color="black", linewidth=.8, sphere=False,
           boundary_color=None, boundary_ls=None, boundary_lw=None,
           mars_lw=0, mars_ls="-"):
    """Draw conics in XY/XZ or surfaces in 3D; YZ requires an X slice in Rm.

    Assumes +X sunward and rotational symmetry about X. YZ shows a section,
    not a projection of the infinite tail. No YZ boundary when x_slice=None.
    """
    plane = plane.upper()
    if plane not in ("XY", "XZ", "YZ", "3D"):
        raise ValueError("plane must be XY, XZ, YZ or 3D")
    color = boundary_color or color
    linewidth = linewidth if boundary_lw is None else boundary_lw
    if sphere:
        plot_mars(ax)
    artists = []
    for name, enabled, ls in (("BS",draw_bs,"--"),("MPB",draw_mpb,":")):
        ls = boundary_ls or ls
        if not enabled:
            continue
        if plane == "YZ":
            if x_slice is None:
                continue
            if not np.isfinite(x_slice):
                raise ValueError("x_slice must be finite")
            k = 0 if name == "BS" or x_slice >= 0 else 1
            x0, length, ecc = MODELS[name][k]
            # r + ecc*(x-x0) = length, r^2=(x-x0)^2+rho^2.
            r = length-ecc*(x_slice-x0)
            rho2 = r*r-(x_slice-x0)**2
            if r > 0 and rho2 >= 0:
                artists.append(ax.add_patch(Circle((0,0),np.sqrt(rho2),fill=False,
                    edgecolor=color,ls=ls,lw=linewidth,label=name)))
            continue
        for x, rho in boundary_profile(name, xmin):
            if plane == "3D":
                phi = np.linspace(0, 2*np.pi, 45)
                xx = np.broadcast_to(x[:,None],(len(x),len(phi)))
                artists.append(ax.plot_wireframe(xx,rho[:,None]*np.cos(phi),
                    rho[:,None]*np.sin(phi),rcount=16,ccount=12,color=color,
                    linewidth=linewidth*.5,alpha=.25))
            else:
                for sign in (1,-1):
                    artists.extend(ax.plot(x,sign*rho,color=color,ls=ls,lw=linewidth))
    return artists


def plot_mars_context(ax, plane="XZ", boundaries=True, **kwargs):
    """Common background used by all trajectory projections."""
    if boundaries:
        bs_mpb(ax,plane=plane,**kwargs)
    plot_mars(ax)
    return ax
