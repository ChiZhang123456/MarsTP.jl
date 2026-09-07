"""Plot/export PSD computed by MarsTP; this script does not bin trajectories.

Usage: python plot_library_psd.py LIBRARY_OUTPUT_ROOT
Input is produced by analyze_saved_probes.jl. PNG and portable sparse NPZ only.
"""
from pathlib import Path
import argparse
import csv
import json
import tomllib
import h5py
import numpy as np
from analyze_monte_carlo import draw_panels


def plot_one(folder):
    meta=tomllib.loads((folder/'library_summary.toml').read_text())
    if (folder/'probe_psd_sparse.npz').exists():
        raise FileExistsError('Derived output exists; choose a fresh analysis run')
    with h5py.File(folder/'library_psd.jld2','r') as f:
        assert bool(f['complete'][()])
        names=['indices_xyz','f3d_s3_m6','number_per_nonzero_bin',
               'unique_particles_per_nonzero_bin','effective_particles_per_nonzero_bin',
               'fxy_s2_m5','fxz_s2_m5','vx_edges_ms','vy_edges_ms','vz_edges_ms','shape_xyz','dv_ms']
        data={name:f[name][()] for name in names}
        # h5py reverses Julia matrix dimensions. COO 3×N becomes N×3;
        # velocity projection matrices must be transposed back to (vx,vy/z).
        data['fxy_s2_m5']=data['fxy_s2_m5'].T
        data['fxz_s2_m5']=data['fxz_s2_m5'].T
        ids=f['particle_ids'][()];rates=f['rate_weights_s'][()];residence=f['residence_s'][()]
    edge=data['vx_edges_ms'];dv=float(data['dv_ms']);density=meta['number_density_m3']
    shape=tuple(data['shape_xyz']);idx=data['indices_xyz']
    assert idx.shape==(len(data['f3d_s3_m6']),3)
    reconstructed=np.zeros(shape[:2])
    np.add.at(reconstructed,(idx[:,0],idx[:,1]),data['f3d_s3_m6']*dv)
    np.testing.assert_allclose(reconstructed,data['fxy_s2_m5'],rtol=1e-12,atol=1e-100)
    for value in (data['f3d_s3_m6'].sum()*dv**3,data['fxy_s2_m5'].sum()*dv**2,data['fxz_s2_m5'].sum()*dv**2):
        assert np.isclose(value,density,rtol=1e-11,atol=1e-100)
    summary={k:meta[k] for k in ('detector_Rm','cube_side_Rm','coordinate_system','source_sha256',
        'source_run_directory','number_density_m3','number_density_cm3','unique_probe_particles',
        'probe_effective_particles','residence_segments','velocity_bin_width_kms','velocity_edges_kms',
        'nonzero_3d_bins','grid_shape','plot_limit_kms','volume_averaged_number_flux_vector_m2_s',
        'density_outside_vlim_m3','estimator','interpolation')}
    summary.update(velocity_range_kms={k:[meta['velocity_min_kms'][i],meta['velocity_max_kms'][i]] for i,k in enumerate('xyz')},
        maximum_flight_age_s=meta['max_flight_time_s'],
        psd_definition='f3d=sum Q_i*tau_ib/(V*dv^3), units s^3/m^6',
        projection_definition='fxy=sum_z f3d*dv; fxz=sum_y f3d*dv, units s^2/m^5',
        storage='COO indices_xyz are zero-based; omitted bins have zero sampled occupancy',
        limitations='Finite Monte Carlo sample and 500 s flight-age window; no steady-state convergence claim.')
    overlap=np.maximum(0,np.minimum(edge[1:],meta['plot_limit_kms']*1000)-np.maximum(edge[:-1],-meta['plot_limit_kms']*1000))
    for name in ('xy','xz'):
        summary[f'displayed_density_fraction_{name}']=float(np.sum(data[f'f{name}_s2_m5']*overlap[:,None]*overlap[None,:])/density) if density else 0.
    events=0
    with (folder/'probe_residence.csv').open(newline='') as src, (folder/'probe_crossings.csv').open('w',newline='') as dst:
        writer=csv.writer(dst)
        writer.writerow(['particle_id','rate_weight_s1','face_flux_m2_s','time_s','face','direction','x_m','y_m','z_m','vx_ms','vy_ms','vz_ms'])
        for row in csv.DictReader(src):
            for ep,key,direction in ((0,'entry_face',1),(1,'exit_face',-1)):
                face=int(row[key])
                if face:
                    q=float(row['rate_weight_s1']);events+=1
                    writer.writerow([row['particle_id'],q,q/meta['cube_face_area_m2'],row[f't{ep}_s'],face,direction,
                        *[row[f'{k}{ep}_m'] for k in 'xyz'],*[row[f'v{k}{ep}_ms'] for k in 'xyz']])
    summary['crossing_events']=events
    with (folder/'probe_particle_weights.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['particle_id','rate_weight_s1','residence_s','density_contribution_m3'])
        writer.writerows(zip(ids,rates,residence,rates*residence/meta['cube_volume_m3']))
    np.savez_compressed(folder/'probe_psd_sparse.npz',**data)
    (folder/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
    draw_panels(folder/'probe_psd_projections',(edge,)*3,data['fxy_s2_m5'],data['fxz_s2_m5'],
        r'Reduced PSD (s$^2$ m$^{-5}$)',
        f'Full velocity integration; {dv/1000:g} km/s bins; maximum flight age {meta["max_flight_time_s"]:g} s',
        meta['unique_probe_particles'],meta['probe_effective_particles'],meta)
    print(f'{folder.name}: {density/1e6:.9g} cm^-3; library PSD exported and plotted')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    root=args.output
    if (root/'library_psd.jld2').exists():
        plot_one(root)
    else:
        assert tomllib.loads((root/'analysis_complete.toml').read_text())['complete']
        for name in ('probe_1_0_2','probe_0_0_2','probe_m1p5_0_1'):plot_one(root/name)


if __name__=='__main__':
    main()
