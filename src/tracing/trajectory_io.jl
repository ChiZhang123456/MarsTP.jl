"""
    write_trajectory_batch(path, trajectories; particle_ids, rate_weights_s,
        source_density_weights_m3=nothing, cell_ids=nothing,
        termination_codes=nothing, species="O2+", coordinate_system="unspecified",
        compress=false)

Write one bounded batch of synchronized Cartesian trajectories to JLD2.
Never overwrites an existing path. Each p<ID>/state is Float64 7×N with rows
t,x,y,z,vx,vy,vz (s,m,m/s); rate_weight_s1 is a physical source rate (s^-1).
Compression is optional. A `complete=true` marker is written only after all
records succeed. Call this after each tracing batch, then release that batch;
do not first collect the entire 50 GB ensemble. Failed partial files are kept
for diagnosis and rejected by the reader. Input accepts trajectories with t/u.
"""
function write_trajectory_batch(path,trajectories;particle_ids,rate_weights_s,
        source_density_weights_m3=nothing,cell_ids=nothing,termination_codes=nothing,
        species="O2+",coordinate_system="unspecified",compress=false)
    ispath(path) && throw(ArgumentError("Refusing to overwrite $path"))
    haskey(TP.SpeciesDict,species) || throw(ArgumentError("Unknown species"))
    n=length(trajectories)
    n>0 && length(particle_ids)==n && length(rate_weights_s)==n || throw(ArgumentError("Batch length mismatch"))
    all(i->i isa Integer && i>0,particle_ids) && length(unique(particle_ids))==n ||
        throw(ArgumentError("Particle IDs must be unique positive integers"))
    all(q->isfinite(q) && q>=0,rate_weights_s) || throw(ArgumentError("Invalid particle rates"))
    for values in (source_density_weights_m3,cell_ids,termination_codes)
        values===nothing || length(values)==n || throw(ArgumentError("Optional batch length mismatch"))
    end
    source_density_weights_m3===nothing || all(q->isfinite(q) && q>=0,source_density_weights_m3) ||
        throw(ArgumentError("Invalid source density weights"))
    cell_ids===nothing || all(i->i isa Integer && i>0,cell_ids) || throw(ArgumentError("Invalid cell IDs"))
    # Validate before opening, so bad input does not create an apparently valid batch.
    for wrapped in trajectories
        traj=_psd_trajectory(wrapped,TP.SpeciesDict[species])
        all(>(0),diff(traj.t)) || throw(ArgumentError("Forward times must strictly increase"))
    end
    mkpath(dirname(abspath(path)))
    jldopen(path,"w";compress) do file
        file["format_version"]=1
        file["particle_ids"]=Int.(particle_ids)
        file["species"]=String(species)
        file["coordinate_system"]=String(coordinate_system)
        file["state_units"]=["s","m","m","m","m s^-1","m s^-1","m s^-1"]
        file["weight_unit"]="s^-1"
        for i in 1:n
            traj=_psd_trajectory(trajectories[i],TP.SpeciesDict[species])
            state=Matrix{Float64}(undef,7,length(traj.t))
            state[1,:].=traj.t
            for j in eachindex(traj.t);state[2:7,j].=traj.u[j][1:6];end
            prefix="p$(particle_ids[i])"
            file["$prefix/state"]=state
            file["$prefix/rate_weight_s1"]=Float64(rate_weights_s[i])
            source_density_weights_m3===nothing || (file["$prefix/source_density_weight_m3"]=Float64(source_density_weights_m3[i]))
            cell_ids===nothing || (file["$prefix/cell_id"]=Int(cell_ids[i]))
            code=termination_codes===nothing ? (hasproperty(traj,:retcode) ? string(traj.retcode) : "unavailable") : string(termination_codes[i])
            file["$prefix/termination_code"]=code
        end
        file["complete"]=true
    end
    return abspath(path)
end

"""
    foreach_saved_trajectory(callback, directory_or_files; species="O2+", progress=nothing)

Read one trajectory at a time and call callback(record), where record contains
particle_id, rate_weight_s, trajectory (t/u), source_density_weight_m3, cell_id,
and termination_code. Supports both write_trajectory_batch format and the
original p<ID>/state legacy batches. A directory must have completion.toml with
complete=true and metadata.toml with matching n_particles/species. Explicit
file lists support partial-run inspection; completeness is then per new batch.
Legacy individual files have no completion marker, reported as legacy format.
The callback must consume the state immediately to keep memory bounded.
"""
function foreach_saved_trajectory(callback,source;species="O2+",progress=nothing)
    haskey(TP.SpeciesDict,species) || throw(ArgumentError("Unknown species"))
    expected=nothing
    if source isa AbstractString && isdir(source)
        completion=TOML.parsefile(joinpath(source,"completion.toml"))
        get(completion,"complete",false)===true || throw(ArgumentError("Incomplete trajectory run"))
        meta=TOML.parsefile(joinpath(source,"metadata.toml"))
        meta["species"]==species || throw(ArgumentError("Saved species mismatch"))
        expected=Int(meta["n_particles"])
        files=sort(filter(p->occursin(r"^trajectories_\d+\.jld2$",basename(p)),readdir(source;join=true)))
    else
        files=source isa AbstractString ? [source] : collect(source)
    end
    isempty(files) && throw(ArgumentError("No trajectory batches"))
    seen=Set{Int}();count=0
    for (batch,path) in enumerate(files)
        jldopen(path,"r") do file
            modern=haskey(file,"format_version")
            if modern
                file["format_version"]==1 || throw(ArgumentError("Unsupported trajectory format"))
                haskey(file,"complete") && file["complete"]===true || throw(ArgumentError("Incomplete batch: $path"))
                file["species"]==species || throw(ArgumentError("Saved species mismatch"))
                file["weight_unit"]=="s^-1" || throw(ArgumentError("Saved rate unit mismatch"))
            end
            ids=modern ? file["particle_ids"] : sort([parse(Int,k[2:end]) for k in keys(file) if occursin(r"^p\d+$",k)])
            for pid in ids
                pid>0 && !(pid in seen) || throw(ArgumentError("Invalid or duplicate saved particle ID $pid"))
                push!(seen,pid)
                prefix="p$pid";state=file["$prefix/state"]
                size(state,1)==7 || throw(ArgumentError("Expected 7×N saved state"))
                ratekey=haskey(file,"$prefix/rate_weight_s1") ? "rate_weight_s1" : "weight_s1"
                rate=Float64(file["$prefix/$ratekey"])
                isfinite(rate) && rate>=0 || throw(ArgumentError("Invalid saved rate"))
                trajectory=(;t=view(state,1,:),u=[view(state,2:7,j) for j in axes(state,2)])
                _psd_trajectory(trajectory,TP.SpeciesDict[species])
                all(>(0),diff(trajectory.t)) || throw(ArgumentError("Invalid saved time sequence"))
                code=haskey(file,"$prefix/termination_code") ? file["$prefix/termination_code"] : "unavailable"
                code in ("unavailable","Success","Terminated","inner","outer","time_limit","zero_rate") ||
                    throw(ArgumentError("Failed/unknown saved termination code: $code"))
                callback((;particle_id=Int(pid),rate_weight_s=rate,trajectory,
                    source_density_weight_m3=haskey(file,"$prefix/source_density_weight_m3") ? file["$prefix/source_density_weight_m3"] : nothing,
                    cell_id=haskey(file,"$prefix/cell_id") ? file["$prefix/cell_id"] : nothing,
                    termination_code=code,legacy_format=!modern))
                count+=1
            end
        end
        progress===nothing || progress((;batch,batches=length(files),particles=count,path))
    end
    expected===nothing || count==expected || throw(ArgumentError("Saved particle count mismatch: $count != $expected"))
    return count
end

"""
    forward_psd_saved(directory_or_files; detector_m, side_m, vlim, vgrid,
        species="O2+", option="3D", velocity_unit=:m_s,
        coordinate_system="unspecified", storage=:dense, progress=nothing)

Streaming equivalent of forward_psd. Reads source rates directly from each
saved trajectory and uses the same accumulator, cube clipping, interpolation,
velocity-bin splitting, and SI normalization. No MHD input or tracing required.
"""
function forward_psd_saved(source;detector_m,side_m,vlim,vgrid,species="O2+",
        option="3D",velocity_unit=:m_s,coordinate_system="unspecified",storage=:dense,progress=nothing)
    acc=ForwardPSDAccumulator(;detector_m,side_m,vlim,vgrid,species,velocity_unit,coordinate_system)
    foreach_saved_trajectory(source;species,progress) do record
        accumulate_forward_psd!(acc,record.trajectory;rate_weight_s=record.rate_weight_s,particle_id=record.particle_id)
        acc.retcodes[end]=record.termination_code
    end
    finish_forward_psd(acc;option,storage)
end
