using MarsTP, StaticArrays, LinearAlgebra

# Reuse the project's spherical-grid moment interpolation; Ui is Cartesian.
mom = MarsTP.load_mhd_moments()
altitudes = isempty(ARGS) ? [200.,400.,600.] : parse.(Float64, ARGS)
all(h -> 200 <= h <= 800, altitudes) || error("Altitude must be in [200,800] km")
scalar(A) = MarsTP.TP.build_interpolator(MarsTP.TP.StructuredGrid, A, mom.r, mom.theta, mom.phi)
println("species,altitude_km,longitude_deg,latitude_deg,n_m3,ux_ms,uy_ms,uz_ms,Ti_K,ur_ms,outward_flux_cm2_s,speed_flux_cm2_s")
for (name, m) in (("O2+", mom.O2plus), ("O+", mom.Oplus))
    source = MarsTP.IonosphereSource(scalar(m.n), scalar(m.t),
        ntuple(k -> scalar(Array(m.v[k,:,:,:])), 3), Rm + 500e3)
    for h in altitudes
      shell = MarsTP.IonosphereSource(source.n, source.temperature, source.velocity, Rm+h*1e3)
      for lat in -90.:2.:90., lon in -180.:2.:180.
        er = SA[cosd(lat)*cosd(lon), cosd(lat)*sind(lon), sind(lat)]
        p = ionosphere_properties(shell, er)
        ur = dot(p.Ui, er)
        println(join((name,h,lon,lat,p.n,p.Ui...,p.Ti,ur,p.n*max(ur,0)/1e4,p.flux/1e4), ','))
      end
    end
end
