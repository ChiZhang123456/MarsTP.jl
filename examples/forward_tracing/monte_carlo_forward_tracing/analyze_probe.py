"""Python regression reference. Production CLI plots MarsTP library output only."""
from pathlib import Path
import argparse
import csv
import json
import tomllib
import numpy as np
from analyze_monte_carlo import load_csv, velocity_segments, draw_panels


def sparse_psd(records, edges, volume):
    nbin = len(edges)-1
    dv = float(edges[1]-edges[0])
    assert np.allclose(np.diff(edges), dv)
    occupancy, particle_bins, particle_total = {}, {}, {}
    current = np.zeros(3)
    for r in records:
        dt = r['t1_s']-r['t0_s']
        assert dt > 0 and r['weight_s1'] >= 0
        pid = int(r['particle_id'])
        amount = r['weight_s1']*dt
        particle_total[pid] = particle_total.get(pid, 0.)+amount
        v0 = np.array([r[f'v{k}0_ms'] for k in 'xyz'])
        v1 = np.array([r[f'v{k}1_ms'] for k in 'xyz'])
        current += amount*(v0+v1)/2/volume
        for idx, fraction in velocity_segments(v0, v1, (edges,)*3):
            key = (idx[0]*nbin+idx[1])*nbin+idx[2]
            value = amount*fraction
            occupancy[key] = occupancy.get(key, 0.)+value
            pk = (pid, key)
            particle_bins[pk] = particle_bins.get(pk, 0.)+value
    keys = np.array(sorted(occupancy), dtype=np.int64)
    values = np.array([occupancy[k] for k in keys])
    indices = np.column_stack((keys//nbin**2, keys//nbin % nbin, keys % nbin)).astype(np.int32)
    psd = values/(volume*dv**3)
    sumsquared, unique = {}, {}
    for (_, key), value in particle_bins.items():
        sumsquared[key] = sumsquared.get(key, 0.)+value**2
        unique[key] = unique.get(key, 0)+1
    neff = np.array([occupancy[k]**2/sumsquared[k] for k in keys])
    n_unique = np.array([unique[k] for k in keys], dtype=np.int32)
    fxy, fxz = np.zeros((nbin,nbin)), np.zeros((nbin,nbin))
    # Sum the stored 3D PSD along the entire omitted axis, then multiply by dv.
    np.add.at(fxy, (indices[:,0],indices[:,1]), psd*dv)
    np.add.at(fxz, (indices[:,0],indices[:,2]), psd*dv)
    density = sum(particle_total.values())/volume
    for recovered in (np.sum(psd)*dv**3, fxy.sum()*dv**2, fxz.sum()*dv**2):
        assert np.isclose(recovered,density,rtol=1e-11,atol=1e-100)
    return dict(indices_xyz=indices, f3d_s3_m6=psd, number_per_nonzero_bin=values,
                unique_particles_per_nonzero_bin=n_unique,
                effective_particles_per_nonzero_bin=neff, fxy_s2_m5=fxy, fxz_s2_m5=fxz), particle_total, current


def analyze_reference(folder, output, dv_kms=5., vmax_kms=500., plot_limit_kms=300.):
    if dv_kms <= 0 or vmax_kms <= 0:
        raise ValueError('Positive grid width and extent required')
    bins = 2*vmax_kms/dv_kms
    if not np.isclose(bins, round(bins)):
        raise ValueError('Velocity span must be an integer multiple of bin width')
    meta = (json.loads((folder/'metadata.json').read_text()) if (folder/'metadata.json').exists()
            else tomllib.loads((folder/'metadata.toml').read_text()))
    records = load_csv(folder/'probe_residence.csv')
    edge = np.linspace(-vmax_kms*1000,vmax_kms*1000,int(round(bins))+1)
    data, contributions, current = sparse_psd(records, edge, meta['cube_volume_m3'])
    total = np.array(list(contributions.values()))
    density = total.sum()/meta['cube_volume_m3']
    neff = total.sum()**2/np.sum(total**2) if len(total) else 0.
    summary = dict(detector_Rm=meta['detector_Rm'], cube_side_Rm=meta['cube_side_Rm'],
        source_run_directory=meta.get('source_run_directory',str(folder.resolve())),
        coordinate_system=meta['coordinate_system'], source_sha256=meta['source_sha256'],
        number_density_m3=float(density), number_density_cm3=float(density/1e6),
        unique_probe_particles=len(contributions), probe_effective_particles=float(neff),
        residence_segments=len(records), velocity_bin_width_kms=dv_kms,
        velocity_edges_kms=[-vmax_kms,vmax_kms], grid_shape=[int(round(bins))]*3,
        nonzero_3d_bins=len(data['f3d_s3_m6']), plot_limit_kms=plot_limit_kms,
        velocity_range_kms={k:[float(min(records[f'v{k}0_ms'].min(),records[f'v{k}1_ms'].min())/1000),
                               float(max(records[f'v{k}0_ms'].max(),records[f'v{k}1_ms'].max())/1000)] for k in 'xyz'} if len(records) else {},
        volume_averaged_number_flux_vector_m2_s=current.tolist(),
        psd_definition='f3d=sum Q_i*tau_ib/(V*dv^3), units s^3/m^6',
        projection_definition='fxy=sum_z f3d*dv; fxz=sum_y f3d*dv, units s^2/m^5',
        storage='COO indices_xyz are zero-based; omitted bins are exactly zero in this Monte Carlo estimator',
        maximum_flight_age_s=meta['max_flight_time_s'],
        limitations='Finite Monte Carlo sample and 500 s flight-age window; no steady-state convergence claim.')
    overlap = np.maximum(0,np.minimum(edge[1:],plot_limit_kms*1000)-np.maximum(edge[:-1],-plot_limit_kms*1000))
    for name in ('xy','xz'):
        summary[f'displayed_density_fraction_{name}'] = float(np.sum(data[f'f{name}_s2_m5']*overlap[:,None]*overlap[None,:])/density) if density else 0.
    output.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(output/'probe_psd_sparse.npz', **data, vx_edges_ms=edge,vy_edges_ms=edge,vz_edges_ms=edge,
                        shape_xyz=np.array(summary['grid_shape']), dv_ms=dv_kms*1000)
    with (output/'probe_particle_weights.csv').open('w',newline='') as f:
        writer=csv.writer(f); writer.writerow(['particle_id','residence_number','density_contribution_m3'])
        for pid, value in sorted(contributions.items()):writer.writerow([pid,value,value/meta['cube_volume_m3']])
    if 'entry_face' in records.dtype.names:
        events = []
        for r in records:
            for endpoint,field,direction in ((0,'entry_face',1),(1,'exit_face',-1)):
                face = int(r[field])
                if face:
                    events.append([int(r['particle_id']),r['weight_s1'],r['weight_s1']/meta['cube_face_area_m2'],
                        r[f't{endpoint}_s'],face,direction,*[r[f'{k}{endpoint}_m'] for k in 'xyz'],
                        *[r[f'v{k}{endpoint}_ms'] for k in 'xyz']])
        with (output/'probe_crossings.csv').open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['particle_id','rate_weight_s1','face_flux_m2_s','time_s','face','direction','x_m','y_m','z_m','vx_ms','vy_ms','vz_ms']);writer.writerows(events)
        summary['crossing_events'] = len(events)
    (output/'analysis_summary.json').write_text(json.dumps(summary,indent=2))
    meta['plot_limit_kms']=plot_limit_kms
    draw_panels(output/'probe_psd_projections',(edge,)*3,data['fxy_s2_m5'],data['fxz_s2_m5'],
        r'Reduced PSD (s$^2$ m$^{-5}$)',
        f'Full velocity integration; {dv_kms:g} km/s bins; maximum flight age {meta["max_flight_time_s"]:g} s',
        len(contributions),neff,meta)
    print(json.dumps(summary,indent=2))
    return summary


def main():
    from plot_library_psd import main as plot_main
    plot_main()


if __name__=='__main__':
    main()
