using StaticArrays, LinearAlgebra
@testset "Omnidirectional DEF from residence and velocity nodes" begin
    m=MarsTP.TP.SpeciesDict["O2+"].m
    alpha=m/(2MarsTP.TP.eV)
    speed=1000.;E=alpha*speed^2
    edges=[0.,E/2,2E,10E]
    # Both endpoints outside: duration in the side-2 cube is 2/speed.
    tr=(;t=[0.,4/speed],u=[SA[-2.,0.,0.,speed,0.,0.],SA[2.,0.,0.,speed,0.,0.]])
    kwargs=(;detector_m=(0.,0.,0.),side_m=2.,rate_weights_s=[3.],energy_edges_eV=edges)
    a=detector_omni_def(tr;kwargs...)
    n=3*(2/speed)/8
    @test a.density_total_m3 ≈ n
    @test a.def[2] ≈ n*E*speed/(4pi*(edges[3]-edges[2]))
    @test a.def[[1,3]]==[0.,0.]
    @test a.units=="eV/(m^2 s eV sr)"
    @test sum(a.dn_dE_m3_eV.*diff(edges))≈n
    # Strongly anisotropic mono-directional VDF and residence agree, without
    # any isotropic assumption. Explicit cell volume is 2*3*4 (m/s)^3.
    b=detector_omni_def(fill(n/24,1,1,1),([speed],[0.],[0.]);
        velocity_cell_widths_m_s=(2.,3.,4.),energy_edges_eV=edges)
    @test b.def≈a.def
    @test b.density_total_m3≈a.density_total_m3
    # DEF ignores the unrelated Cartesian velocity limits and projection.
    p=forward_psd(tr;kwargs...,vlim=1.,vgrid=2,option="Vx-Vz")
    @test p.omni_def.def≈a.def
    @test p.density_in_range_m3==0
    @test p.density_total_m3≈n
    @test forward_psd(tr;detector_m=(0.,0.,0.),side_m=2.,rate_weights_s=[3.],vlim=2000.,vgrid=2).omni_def===nothing
    # Reversal crosses each positive energy shell twice. Integral |v|^3 over
    # [-speed,speed] gives exact Ev-weighted interval contributions.
    turn=(;t=[0.,2.],u=[SA[0.,0.,0.,-speed,0.,0.],SA[0.,0.,0.,speed,0.,0.]])
    d=detector_omni_def(turn;detector_m=(0.,0.,0.),side_m=2.,rate_weights_s=[4.],energy_edges_eV=[0.,E/4,E])
    expected1=(4/8)*alpha*speed^3/32
    expected2=(4/8)*alpha*speed^3*15/32
    @test d.def.*(4pi.*diff(d.energy_edges_eV))≈[expected1,expected2] rtol=1e-9
    @test d.density_per_bin_m3≈[.5,.5]
    # Outside-energy density is reported and not renormalized into valid bins.
    out=detector_omni_def(tr;kwargs...,energy_edges_eV=[0.,E/2])
    @test out.def==[0.]
    @test out.density_outside_energy_range_m3≈n
    # Final upper edge is included; zero speed contributes density but no DEF.
    lastbin=detector_omni_def(tr;kwargs...,energy_edges_eV=[0.,E])
    @test lastbin.density_in_range_m3≈n
    still=(;t=[0.,1.],u=[zeros(6),zeros(6)])
    z=detector_omni_def(still;kwargs...)
    @test z.def==zeros(3)
    @test z.density_per_bin_m3[1]≈3/8
    # Smooth non-collinear speed variation, checked against high-resolution
    # composite trapezoid integration independently of the adaptive quadrature.
    changing=(;t=[0.,1.],u=[SA[0.,0.,0.,-1000.,200.,0.],SA[0.,0.,0.,1500.,-50.,400.]])
    smooth=detector_omni_def(changing;kwargs...,energy_edges_eV=[0.,100.])
    vv(t)=SA[-1000.,200.,0.]+t*SA[2500.,-250.,400.]
    fun(t)=alpha*norm(vv(t))^3
    steps=20000
    reference=(sum(fun(i/steps) for i in 1:steps-1)+(fun(0.)+fun(1.))/2)/steps
    @test smooth.def[1]≈3/8*reference/(4pi*100.) rtol=2e-8
    # Shared backward slice accumulation equals full 3D rectangular quadrature.
    axes=([-speed,0.,speed],[-speed,speed],[0.,speed])
    f=reshape(collect(1.:12.),3,2,2)
    e=[0.,E/2,1.5E,4E]
    full=detector_omni_def(f,axes;velocity_cell_widths_m_s=(speed,2speed,speed),energy_edges_eV=e)
    acc=MarsTP._OmniDEFAccumulator(e,"O2+")
    for iy in 1:2
        MarsTP._def_vdf!(acc,reshape(f[:,iy,:],3,1,2),(axes[1],[axes[2][iy]],axes[3]),(speed,2speed,speed))
    end
    @test detector_omni_def(acc).def≈full.def
    @test detector_omni_def(acc).density_per_bin_m3≈full.density_per_bin_m3
    # Saved-trajectory API uses the same DEF accumulator.
    path=tempname()*".jld2"
    write_trajectory_batch(path,[tr];particle_ids=[1],rate_weights_s=[3.],species="O2+")
    saved=forward_psd_saved([path];detector_m=(0.,0.,0.),side_m=2.,vlim=1.,vgrid=2,energy_edges_eV=edges)
    @test saved.omni_def.def≈a.def
    @test_throws ArgumentError detector_omni_def(tr;kwargs...,energy_edges_eV=[1.,1.])
    @test_throws ArgumentError detector_omni_def(tr;kwargs...,energy_edges_eV=[-1.,1.])
    @test_throws ArgumentError detector_omni_def(tr;kwargs...,rate_weights_s=[-1.])
    @test_throws ArgumentError detector_omni_def((;t=reverse(tr.t),u=reverse(tr.u));kwargs...)
    @test_throws ArgumentError detector_omni_def(-ones(1,1,1),([speed],[0.],[0.]);velocity_cell_widths_m_s=(1.,1.,1.),energy_edges_eV=edges)
end
