"""Analytic SI normalization and correlated-sample tests; no MHD data needed."""
import unittest
import numpy as np
from analyze_monte_carlo import bin_residence, velocity_segments


def record(pid, weight, t0, t1, v0, v1):
    return dict(particle_id=pid, weight_s1=weight, t0_s=t0, t1_s=t1,
                **{f"v{k}0_ms": v0[i] for i,k in enumerate("xyz")},
                **{f"v{k}1_ms": v1[i] for i,k in enumerate("xyz")})


class TestPSD(unittest.TestCase):
    def test_bin_crossing(self):
        edges=(np.array([-2.,0.,2.]),)*3
        pieces=list(velocity_segments(np.array([-1.,1.,1.]),np.array([1.,1.,1.]),edges))
        self.assertEqual(pieces,[((0,1,1),.5),((1,1,1),.5)])

    def test_density_units_and_sample_correlation(self):
        edges=(np.array([0.,10.]),)*3
        v=np.array([2.,3.,4.])
        rows=[record(1,100.,0.,1.,v,v),record(1,100.,1.,2.,v,v)]
        f,number,unique,neff,j=bin_residence(rows,edges,volume=20.)
        self.assertAlmostEqual(number.item(),200.)
        self.assertAlmostEqual(f.item(),.01)  # 200/(20*10^3), SI PSD
        self.assertEqual(unique.item(),1)
        self.assertEqual(neff.item(),1.)
        np.testing.assert_allclose(j,10*v)

    def test_independent_weights(self):
        edges=(np.array([0.,10.]),)*3
        v=np.array([2.,3.,4.])
        rows=[record(1,1.,0.,1.,v,v),record(2,3.,0.,1.,v,v)]
        _,_,unique,neff,_=bin_residence(rows,edges,volume=1.)
        self.assertEqual(unique.item(),2)
        self.assertAlmostEqual(neff.item(),1.6)

    def test_outside_bin_is_not_silently_dropped(self):
        edges=(np.array([0.,10.]),)*3
        with self.assertRaises(ValueError):
            list(velocity_segments(np.array([11.,1.,1.]),np.array([12.,1.,1.]),edges))


if __name__=="__main__":
    unittest.main()
