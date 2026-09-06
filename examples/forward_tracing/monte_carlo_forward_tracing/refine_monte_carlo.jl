# Deterministic subset validation, with identical initial states and weights.
# Usage: julia --project=. refine_monte_carlo.jl RUN_DIR
include("monte_carlo_shell.jl")
using .ShellMonteCarlo, MarsTP, StaticArrays, LinearAlgebra, TOML
const MC=ShellMonteCarlo
function main(out)
    isfile(joinpath(out,"completion.toml")) || error("Full run must complete first")
    meta=TOML.parsefile(joinpath(out,"metadata.toml"))
    if !haskey(meta,"particle_mass_kg")
        meta["particle_mass_kg"]=MarsTP.TP.SpeciesDict["O2+"].m
        meta["particle_charge_C"]=MarsTP.TP.SpeciesDict["O2+"].q
        open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,meta);end
    end
    names=split(readline(joinpath(out,"particles.csv")),',');lookup=Dict(String(n)=>i for (i,n) in enumerate(names))
    haskey(lookup,"rate_weight_s1") && (lookup["weight_s1"]=lookup["rate_weight_s1"])
    val(row,name)=parse(Float64,row[lookup[name]])
    spaced=Set(round.(Int,range(1,meta["n_particles"],length=64)))
    energy_mode=length(ARGS)>1 && ARGS[2]=="energy"
    data=Vector{SubString{String}}[]
    scores=Float64[]
    all_hits=0
    open(joinpath(out,"particles.csv")) do io
        readline(io)
        for (i,line) in enumerate(eachline(io))
            r=split(line,',')
            ishit=val(r,"probe_residence_s")>0
            all_hits+=ishit
            val(r,"weight_s1")>0 || continue
            if energy_mode
                score=abs(val(r,"residual_eV"))/max(abs(val(r,"deltaK_eV")),abs(val(r,"work_eV")),1.)
                if length(scores)<16
                    push!(data,r);push!(scores,score)
                elseif score>minimum(scores)
                    j=argmin(scores);data[j]=r;scores[j]=score
                end
            elseif ishit || i in spaced
                push!(data,r)
            end
        end
    end
    # All probe contributors (up to 200 highest density contributions), plus
    # 32 spatially spaced draws. Selected probe subset is not the full source.
    hit=findall(r->val(r,"probe_residence_s")>0,data)
    sort!(hit;by=i->val(data[i],"probe_residence_s")*val(data[i],"weight_s1"),rev=true)
    selected=energy_mode ? collect(eachindex(data)) : sort!(unique(vcat(hit[1:min(200,length(hit))],findall(r->val(r,"probe_residence_s")==0,data))))
    fields=load_mhd_fields()
    c=MC.Config()
    fields,source,_=MC.validate_geometry(fields.path,fields,c)
    param=MarsTP.mhd_param(fields;species="O2+")
    target=joinpath(out,length(ARGS)>1 && ARGS[2]=="energy" ? "energy_refinement.csv" : "timestep_refinement.csv")
    ispath(target) && error("Refinement output already exists")
    open(target,"w") do io
        println(io,"particle_id,weight_s1,dt_s,status,end_time_s,x_m,y_m,z_m,vx_ms,vy_ms,vz_ms,probe_residence_s,probe_crossings,energy_residual_eV,deltaK_eV")
        for i in selected
            r=data[i]
            id=parse(Int,r[lookup["particle_id"]]);W=val(r,"weight_s1")
            x=SA[val(r,"x0_m"),val(r,"y0_m"),val(r,"z0_m")]
            v=SA[val(r,"vx0_ms"),val(r,"vy0_ms"),val(r,"vz0_ms")]
            for dt in (.1,.05,.025)
                result=MC.trace_particle(x,v,param,MC.Config(dt=dt,tmax=meta["max_flight_time_s"]);keep_history=false)
                dwell=sum((a[2]-a[1] for a in result.residence);init=0.)
                if dt==.1
                    @assert result.status==r[lookup["status"]]
                    @assert isapprox(dwell,val(r,"probe_residence_s");rtol=1e-10,atol=1e-10)
                end
                println(io,join((id,W,dt,result.status,result.time,result.x...,result.v...,dwell,length(result.events),result.residual,result.dK),','))
            end
        end
    end
    open(joinpath(out,energy_mode ? "energy_refinement_selection.toml" : "refinement_selection.toml"),"w") do io
        TOML.print(io,Dict("particles_refined"=>length(selected),"all_coarse_probe_contributors"=>all_hits,
            "selection"=>energy_mode ? "16 largest relative energy residuals" : "up to 200 largest probe density contributors plus positive-rate spaced controls"))
    end
    println("Refined $(length(selected)) deterministic particles; total coarse probe contributors=$all_hits.")
end
main(ARGS[1])
