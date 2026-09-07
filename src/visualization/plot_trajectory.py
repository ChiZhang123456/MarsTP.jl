"""Plot forward/backward trajectory files in 2D or 3D, without py_space_zc."""
from pathlib import Path
import argparse
import json
import tomllib
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
try:
    from .mars import plot_mars_context
except ImportError:
    from mars import plot_mars_context


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def _metadata(folder):
    for name in ("metadata.toml", "metadata.json"):
        path = folder/name
        if path.exists():
            raw = path.read_text(encoding="utf8")
            return tomllib.loads(raw) if path.suffix == ".toml" else json.loads(raw)
    return {}


def load_trajectories(source, *, species=None, assume_species=None, count=5000,
                      max_points=300, seed=20260906, position_unit=None, rm_m=3390000.,
                      particle_ids=None):
    """Reservoir sample trajectories, returning points in Rm.

    Supports MarsTP p<ID>/state JLD2/HDF5 (N x 7 through h5py), and JSONL
    records with id, points (N x 3), times_s. JSONL position units must be
    declared using position_unit or a record/sidecar position_unit field.
    Species filtering uses record, file or sidecar metadata. assume_species
    explicitly labels old data without a species field; never overrides one.
    Handles increasing forward or decreasing backward times, preserving order.
    """
    if count < 1 or max_points < 2 or not np.isfinite(rm_m) or rm_m <= 0:
        raise ValueError("count >= 1, max_points >= 2, rm_m > 0 required")
    wanted = None if particle_ids is None else set(particle_ids)
    source = Path(source)
    if source.is_dir():
        files = sorted(source.glob("trajectories_*.jld2"))
        if not files and (source/"trajectories.jsonl").exists():
            files = [source/"trajectories.jsonl"]
    else:
        files = [source]
    if not files:
        raise ValueError("No trajectory files found (PSD-only files contain no paths)")
    rng, selected, eligible = np.random.default_rng(seed), [], 0

    def consider(item):
        nonlocal eligible
        if wanted is not None and item["id"] not in wanted:
            return
        known = item["species"] or assume_species
        if species is not None and known is None:
            raise ValueError("Species metadata missing; supply assume_species explicitly")
        if species is not None and known != species:
            return
        item["species"] = known or "unspecified"
        eligible += 1
        if len(selected) < count:
            selected.append(item)
        else:
            j = int(rng.integers(eligible))
            if j < count:
                selected[j] = item

    for path in files:
        meta = _metadata(path.parent)
        if path.suffix in (".jld2", ".h5", ".hdf5"):
            if position_unit not in (None, "m"):
                raise ValueError("MarsTP state format has fixed SI position units (m)")
            with h5py.File(path,"r") as f:
                if "format_version" in f and ("complete" not in f or not bool(f["complete"][()])):
                    raise ValueError(f"Incomplete trajectory batch: {path}")
                kind = _text(f["species"][()]) if "species" in f else meta.get("species")
                frame = _text(f["coordinate_system"][()]) if "coordinate_system" in f else meta.get("coordinate_system","unspecified")
                found = False
                for key in f:
                    if not (key.startswith("p") and key[1:].isdigit() and isinstance(f[key],h5py.Group) and "state" in f[key]):
                        continue
                    found = True
                    state = f[key+"/state"]
                    if state.ndim != 2 or state.shape[1] != 7:
                        raise ValueError(f"Expected N x 7 SI state: {path}:{key}")
                    if state.shape[0] < 2:
                        continue
                    record_species = _text(f[key+"/species"][()]) if key+"/species" in f else kind
                    consider(dict(path=path,key=key,id=int(key[1:]),species=record_species,
                                  coordinate_system=frame,unit="m"))
                if not found:
                    raise ValueError(f"No saved trajectories in {path}; PSD arrays cannot reconstruct paths")
        elif path.suffix == ".jsonl":
            with path.open(encoding="utf8") as stream:
                for line in stream:
                    r = json.loads(line)
                    unit = position_unit or r.get("position_unit") or meta.get("position_unit")
                    if unit not in ("m","km","Rm"):
                        raise ValueError("JSONL needs explicit position_unit: m, km or Rm")
                    # Bound reservoir memory even for very long JSONL paths.
                    points, times = np.asarray(r["points"],float), np.asarray(r["times_s"],float)
                    if points.shape != (len(times),3) or len(times)<2:
                        raise ValueError("Expected N x 3 points and matching times_s")
                    ix = np.unique(np.linspace(0,len(times)-1,min(max_points,len(times)),dtype=int))
                    r["points"], r["times_s"] = points[ix], times[ix]
                    consider(dict(path=path,id=r["id"],species=r.get("species",meta.get("species")),
                        coordinate_system=r.get("coordinate_system",meta.get("coordinate_system","unspecified")),
                        unit=unit,raw=r))
        else:
            raise ValueError(f"Unsupported trajectory format: {path.suffix}")
    if not selected:
        raise ValueError("No trajectories match the requested species")
    result = []
    # Only selected HDF5 trajectories are read; bounded display memory.
    for item in selected:
        if "key" in item:
            with h5py.File(item["path"],"r") as f:
                ds = f[item["key"]+"/state"]
                ix = np.unique(np.linspace(0,len(ds)-1,min(max_points,len(ds)),dtype=int))
                a = ds[ix,:4]
                times, xyz = a[:,0], a[:,1:4]
        else:
            r = item.pop("raw")
            xyz, times = np.asarray(r["points"],float), np.asarray(r["times_s"],float)
            if xyz.ndim != 2 or xyz.shape != (len(times),3) or len(times)<2:
                raise ValueError("Expected matching N x 3 points and times_s, N >= 2")
            ix = np.unique(np.linspace(0,len(times)-1,min(max_points,len(times)),dtype=int))
            xyz,times = xyz[ix],times[ix]
        if not np.isfinite(xyz).all() or not np.isfinite(times).all():
            raise ValueError("Nonfinite trajectory")
        dt = np.diff(times)
        if not (np.all(dt>0) or np.all(dt<0)):
            raise ValueError("Trajectory times must be strictly monotonic")
        scale = {"m":1/rm_m,"km":1000/rm_m,"Rm":1}[item["unit"]]
        result.append(dict(id=item["id"],source=str(item["path"]),points=xyz*scale,
            times_s=times,species=item["species"],coordinate_system=item["coordinate_system"]))
    return result


def plot_trajectory(source, *, planes=("XZ","XY","YZ"), limits=(-4.2,4.2),
                    boundaries=True, x_slice=None, output=None, **load_options):
    """Return (figure, axes, sampled_records). planes may include '3D'.

    BS/MPB assume +X sunward; disable boundaries for other coordinate frames.
    Particle species is metadata selection, not a mass/charge transformation.
    """
    records = load_trajectories(source,**load_options)
    planes = (planes,) if isinstance(planes,str) else tuple(planes)
    planes = tuple(p.upper() for p in planes)
    if not planes or any(p not in ("XY","XZ","YZ","3D") for p in planes):
        raise ValueError("planes: XY, XZ, YZ, 3D")
    if len(limits)!=2 or not np.isfinite(limits).all() or limits[0]>=limits[1]:
        raise ValueError("limits must be finite and increasing")
    frames = {r["coordinate_system"] for r in records}
    if len(frames)>1:
        raise ValueError("Cannot overlay trajectories with different coordinate frames")
    with mpl.rc_context({"font.family":"Arial","font.size":10}):
        fig = plt.figure(figsize=(5*len(planes),5),layout="constrained")
        axes = []
        colors = plt.get_cmap("turbo")(np.linspace(.05,.95,len(records)))
        for k,plane in enumerate(planes):
            ax = fig.add_subplot(1,len(planes),k+1,projection="3d" if plane=="3D" else None)
            axes.append(ax)
            plot_mars_context(ax,plane,boundaries=boundaries,xmin=limits[0],x_slice=x_slice)
            for r,color in zip(records,colors):
                p = r["points"]
                if plane=="3D":
                    ax.plot(*p.T,color=color,lw=.5,alpha=.4)
                else:
                    i,j = ["XYZ".index(v) for v in plane]
                    ax.plot(p[:,i],p[:,j],color=color,lw=.5,alpha=.4,zorder=2)
            if plane=="3D":
                ax.set(xlim=limits,ylim=limits,zlim=limits,xlabel="X (Rm)",ylabel="Y (Rm)",zlabel="Z (Rm)")
                ax.set_box_aspect((1,1,1))
            else:
                ax.set(xlim=limits,ylim=limits,xlabel=plane[0]+" (Rm)",ylabel=plane[1]+" (Rm)",aspect="equal")
            ax.set_title(plane + (f" (boundaries: X={x_slice:g} Rm)" if plane=="YZ" and x_slice is not None else ""))
        fig.suptitle(", ".join(sorted({r["species"] for r in records}))+"; "+next(iter(frames)))
        if output is not None:
            output=Path(output)
            if output.suffix.lower()!=".png":
                raise ValueError("Output must be PNG")
            output.parent.mkdir(parents=True,exist_ok=True)
            fig.savefig(output,dpi=180)
    return fig,axes,records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source")
    parser.add_argument("--output",required=True)
    parser.add_argument("--species")
    parser.add_argument("--particle-ids",nargs="+",type=int)
    parser.add_argument("--assume-species")
    parser.add_argument("--position-unit",choices=("m","km","Rm"))
    parser.add_argument("--count",type=int,default=5000)
    parser.add_argument("--max-points",type=int,default=300)
    parser.add_argument("--seed",type=int,default=20260906)
    parser.add_argument("--rm-m",type=float,default=3390000.)
    parser.add_argument("--planes",nargs="+",default=["XZ","XY","YZ"])
    parser.add_argument("--x-slice",type=float)
    parser.add_argument("--no-boundaries",action="store_true")
    args=vars(parser.parse_args())
    args["boundaries"]=not args.pop("no_boundaries")
    fig,_,records=plot_trajectory(**args)
    plt.close(fig)
    print(f"Saved {len(records)} trajectories")


if __name__=="__main__":
    main()
