# Reconstruct the original Boris within-step synchronized velocities.
# No trajectories are reintegrated. Usage: julia --project=. synchronize_reprobe.jl OUTPUT
include("monte_carlo_shell.jl")
using .ShellMonteCarlo, MarsTP, TOML, SHA
const MC = ShellMonteCarlo
function main(root)
    isfile(joinpath(root,"extraction_complete.json")) || error("Incomplete extraction")
    meta=TOML.parsefile(joinpath(root,"source_metadata.toml"))
    fields = load_mhd_fields()
    fields, _, _ = MC.validate_geometry(fields.path, fields, MC.Config())
    param = MarsTP.mhd_param(fields;species="O2+")
    bytes2hex(open(sha256,fields.path)) == meta["source_sha256"] || error("MHD input differs from the saved run")
    dt = meta["dt_s"]
    for name in ("probe_1_0_2", "probe_0_0_2", "probe_m1p5_0_1")
        folder=joinpath(root,name)
        target=joinpath(folder,"probe_residence.csv")
        isfile(target) && error("Refusing to overwrite $target")
        open(target,"w") do out
            println(out,"particle_id,rate_weight_s1,t0_s,t1_s,x0_m,y0_m,z0_m,x1_m,y1_m,z1_m,vx0_ms,vy0_ms,vz0_ms,vx1_ms,vy1_ms,vz1_ms,entry_face,exit_face")
            for (i,line) in enumerate(eachline(joinpath(folder,"clipped_segments.csv")))
                i==1 && continue
                a=parse.(Float64,split(line,','))
                x0,x1=MC.Vec(a[5:7]),MC.Vec(a[8:10])
                half=MC.Vec(a[12:14]); tbase=a[11]
                0<=a[3]-tbase<=dt+1e-10 || error("Unexpected saved step")
                0<a[4]-tbase<=dt+1e-10 || error("Unexpected saved step")
                v0=MC.TP.update_velocity(half,x0,a[3]-tbase-dt/2,a[3],param)
                v1=MC.TP.update_velocity(half,x1,a[4]-tbase-dt/2,a[4],param)
                all(isfinite,v0) && all(isfinite,v1) || error("Nonfinite velocity")
                println(out,join((Int(a[1]),a[2:10]...,v0...,v1...,Int(a[15]),Int(a[16])),','))
            end
        end
        println("Synchronized $name")
    end
end
main(ARGS[1])
