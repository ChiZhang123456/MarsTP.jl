using MarsTP, StaticArrays, LinearAlgebra
# Sample the exact same MHD moment interpolators used by backtracing.
source=load_ionosphere_source(;altitude_km=200.)
println("altitude_km,longitude_deg,latitude_deg,n_m3,ux_ms,uy_ms,uz_ms,flux_m2_s,Ti_K")
for h in (200.,400.)
    shell=IonosphereSource(source.n,source.temperature,source.velocity,Rm+h*1e3)
    for lat in -90.:2.:90., lon in -180.:2.:180.
        d=SA[cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat)]
        m=ionosphere_properties(shell,d)
        println(join((h,lon,lat,m.n,m.Ui...,m.flux,m.Ti),','))
    end
end
