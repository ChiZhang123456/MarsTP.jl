const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))

project_path(parts...) = joinpath(PROJECT_ROOT, parts...)
data_path(parts...) = project_path("data", parts...)

function resolve_project_path(path::AbstractString)
    isabspath(path) && return path
    isfile(path) && return abspath(path)
    return project_path(path)
end
