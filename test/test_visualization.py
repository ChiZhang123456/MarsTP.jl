"""Run: python test/test_visualization.py (no simulation required)."""
import sys
import json
import tempfile
from pathlib import Path
import unittest
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from visualization import boundary_profile, bs_mpb, load_trajectories, plot_trajectory


class VisualizationTest(unittest.TestCase):
    def test_conics(self):
        from visualization.mars import MODELS
        for name in MODELS:
            for (x,rho),(x0,length,ecc) in zip(boundary_profile(name),MODELS[name]):
                np.testing.assert_allclose(np.sqrt((x-x0)**2+rho**2)+ecc*(x-x0),length,atol=1e-12)
                self.assertTrue(np.all(rho>=0))
        fig,ax=plt.subplots()
        circles=bs_mpb(ax,plane="YZ",x_slice=0)
        self.assertEqual(len(circles),2)
        plt.close(fig)

    def test_read_and_plot(self):
        # Retain tiny fixtures in the output folder; no deletion of user files.
        root=Path(__file__).resolve().parents[1]/"outputs"
        root.mkdir(exist_ok=True)
        folder=Path(tempfile.mkdtemp(prefix="visualization_test_",dir=root))
        path=folder/"trajectories_1.jld2"
        with h5py.File(path,"w") as f:
            f["format_version"]=1; f["complete"]=True; f["species"]="O2+"
            f["coordinate_system"]="MSO"
            a=np.zeros((4,7));a[:,0]=np.arange(4);a[:,1]=np.arange(4)*3390000
            f["p1/state"]=a
            f["p2/state"]=a;f["p2/species"]="O+"
        records=load_trajectories(folder,species="O2+",max_points=2)
        self.assertEqual(len(records),1)
        np.testing.assert_allclose(records[0]["points"][:,0],[0,3])
        self.assertEqual(load_trajectories(folder,species="O+")[0]["id"],2)
        with self.assertRaises(ValueError):load_trajectories(folder,species="H+")
        j=folder/"backward.jsonl"
        j.write_text(json.dumps(dict(id=1,points=[[0,0,2],[1,0,2]],times_s=[0,-1]))+"\n")
        with self.assertRaises(ValueError):load_trajectories(j,position_unit="Rm",species="O2+")
        back=load_trajectories(j,position_unit="Rm",species="O2+",assume_species="O2+")
        self.assertEqual(back[0]["times_s"][-1],-1)
        fig,axes,_=plot_trajectory(path,planes=("XZ","XY","YZ","3D"),x_slice=0,
                                  output=folder/"smoke.png")
        self.assertEqual(len(axes),4);plt.close(fig)
        with h5py.File(path,"r+") as f:f["complete"][...]=False
        with self.assertRaises(ValueError):load_trajectories(path)


if __name__=="__main__":unittest.main()
