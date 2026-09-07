"""Geometry and sparse/dense normalization equivalence, independent of MHD."""
import unittest
import numpy as np
from reprobe_saved import intersections
from analyze_probe import sparse_psd
from analyze_monte_carlo import bin_residence
from test_monte_carlo_analysis import record


class TestReprobe(unittest.TestCase):
    def test_cube_entry_exit(self):
        state=np.array([[0.,-2.,0.,0.,4.,0.,0.],[1.,2.,0.,0.,4.,0.,0.]])
        hit=intersections(state,np.zeros((1,3)),2.)
        self.assertEqual(hit,[(0,0,.25,.75,-1,1)])

    def test_parallel_miss(self):
        state=np.array([[0.,-2.,2.,0.,4.,0.,0.],[1.,2.,2.,0.,4.,0.,0.]])
        self.assertEqual(intersections(state,np.zeros((1,3)),2.),[])

    def test_sparse_dense_and_units(self):
        rows=[record(1,100.,0.,1.,[-1500.,500.,500.],[1500.,500.,500.]),
              record(1,100.,1.,2.,[1500.,500.,500.],[1500.,500.,500.]),
              record(2,50.,0.,1.,[500.,500.,500.],[500.,500.,500.])]
        edge=np.arange(-2000.,2001.,1000.)
        sparse,total,current=sparse_psd(rows,edge,20.)
        dense,number,unique,neff,j=bin_residence(rows,(edge,)*3,20.)
        idx=tuple(sparse['indices_xyz'].T)
        np.testing.assert_allclose(sparse['f3d_s3_m6'],dense[idx],rtol=1e-14)
        np.testing.assert_allclose(sparse['fxy_s2_m5'],dense.sum(axis=2)*1000,rtol=1e-14)
        np.testing.assert_allclose(sparse['fxz_s2_m5'],dense.sum(axis=1)*1000,rtol=1e-14)
        np.testing.assert_allclose(sparse['effective_particles_per_nonzero_bin'],neff[idx])
        np.testing.assert_array_equal(sparse['unique_particles_per_nonzero_bin'],unique[idx])
        np.testing.assert_allclose(current,j)
        self.assertAlmostEqual(sparse['fxy_s2_m5'].sum()*1000**2,12.5)
        self.assertEqual(total,{1:200.,2:50.})

    def test_outside_grid_raises(self):
        with self.assertRaises(ValueError):
            sparse_psd([record(1,1.,0.,1.,[3.,0.,0.],[3.,0.,0.])],np.array([-1.,0.,1.]),1.)


if __name__=='__main__':
    unittest.main()
