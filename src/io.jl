function enforce_surface!(data, input, model, grid; effective_radius=input.sradius)
    traction = zeros(input.nt, length(input.r))
    if input.istype == 0 && input.zstop == 0 && input.zsbtm == 0 && input.isurfcon > 0
        w = [r < effective_radius ? 2(1-(r/effective_radius)^2)/(pi*effective_radius^2) : 0.0 for r in input.r]
        times = range(0.0, input.twindow; length=input.nt)
        traction .= -source_value.(Ref(input), times) .* transpose(w)
    end
    for (iz, l) in enumerate(grid.lzr)
        l == 1 && input.isurfcon != 0 || continue
        la, mu, alf = model.la[1], model.mu[1], model.alf[1]
        data[:, :, iz, 6] .= 0
        if input.isurfcon == 1
            data[:, :, iz, 8] .= 0
        else
            data[:, :, iz, 9] .= 0
        end
        data[:, :, iz, 3] .= (-la .* (data[:, :, iz, 4] + data[:, :, iz, 5]) +
                              alf .* data[:, :, iz, 8] + traction) ./ (la + 2mu)
    end
    data
end

"""Save a complete JLD2 archive by default; use :text or :both for analysis tables."""
function write_outputs(result; output_dir::AbstractString="results/simulation", format=:jld2)
    format = Symbol(format)
    format in (:text, :jld2, :both) || error("Output format must be jld2, text, or both.")
    input, data = result.input, result.data
    mkpath(output_dir)
    nt, nr, nzr, _ = size(data)
    dt = input.twindow/(nt-1)
    paths = String[]
    if format != :text
        push!(paths, save_result(joinpath(output_dir, "result.jld2"), result))
        format == :jld2 && return paths
    end
    labels = ("Uz", "Ur", "Ezz", "Err", "Ett", "Ezr", "Tlt", "Pp", "Dvz", "Dvr")
    for ic in 1:10
        input.selected[ic] || continue
        path = joinpath(output_dir, input.filenames[ic])
        mkpath(dirname(path))
        open(path, "w") do io
            print(io, "Time[s]")
            for iz in 1:nzr, ir in 1:nr
                @printf(io, " %s_Z%04d_R%04d", labels[ic], iz, ir)
            end
            println(io)
            for it in 1:nt
                @printf(io, "%.16e", result.time[it])
                for iz in 1:nzr, ir in 1:nr
                    @printf(io, " %.16e", data[it, ir, iz, ic])
                end
                println(io)
            end
        end
        push!(paths, abspath(path))
    end
    for (time, filename) in zip(input.timesn, input.filesn)
        path = joinpath(output_dir, filename)
        mkpath(dirname(path))
        plane = snapshot(result, time).data
        open(path, "w") do io
            println(io, "Depth[m] Distance[m] Uz Ur Ezz Err Ett Ezr Tlt Pp Dvz Dvr")
            for iz in 1:nzr, ir in 1:nr
                @printf(io, "%.16e %.16e", input.zr[iz], input.r[ir])
                for ic in 1:10
                    @printf(io, " %.16e", plane[ir, iz, ic])
                end
                println(io)
            end
        end
        push!(paths, abspath(path))
    end
    if result.equivalent !== nothing
        path = joinpath(output_dir, input.equivalent_filename)
        mkpath(dirname(path))
        open(path, "w") do io
            println(io, input.istype == 2 ? "Time[s] Mean_source_pressure_volume[Pa*m^3]" : "Time[s] Injection_rate[m^3/s]")
            for it in 1:nt
                @printf(io, "%.16e %.16e\n", result.time[it], result.equivalent[it])
            end
        end
        push!(paths, abspath(path))
    end
    path = joinpath(output_dir, "simulation_info.txt")
    open(path, "w") do io
        println(io, "MaoyePoro ", pkgversion(@__MODULE__), "; coupled poroelastic deformation and fluid diffusion")
        println(io, "SI units; z downward; source type: ", input.istype,
            "; surface condition: ", input.isurfcon)
        println(io, "source meaning: ", SOURCE_KINDS[input.istype+1], " = ", SOURCE_DESCRIPTIONS[input.istype+1])
        if result.equivalent !== nothing
            println(io, input.istype == 2 ?
                "equivalent output: mean source pore pressure * source volume [Pa*m^3 = N*m]" :
                "equivalent output: fluid volume injection rate [m^3/s]; negative values mean extraction; pulse volumes are recorded separately")
        end
        println(io, "source screen [m]: ", (input.zstop, input.zsbtm), "; input radius [m]: ", input.sradius)
        println(io, "receiver depths [m]: ", input.zr, "; radii [m]: ", input.r)
        println(io, "time samples: ", nt, "; dt [s]: ", dt, "; window [s]: ", input.twindow)
        println(io, "homogeneous layers: ", length(result.model.h),
            "; depth levels: ", length(result.grid.z), "; requested integration accuracy: ", input.accuracy)
        for d in result.diagnostics
            println(io, d)
        end
        result.equivalent_diagnostics === nothing || println(io, "equivalent source: ", result.equivalent_diagnostics)
        hasproperty(result, :settings) && println(io, "solver settings: ", result.settings)
        println(io, "Convergence flags concern numerical truncation criteria; they are not an error bound for the physical model.")
    end
    push!(paths, abspath(path))
    paths
end

"""Read field time series; columns are time, then depth-major/radius-minor samples."""
function read_timeseries(path::AbstractString)
    open(path) do io
        names = split(readline(io))
        values = Float64[]
        for line in eachline(io)
            tokens = split(line)
            isempty(tokens) && continue
            length(tokens) == length(names) || error("Time-series row has the wrong number of columns.")
            for token in tokens
                push!(values, parse(Float64, token))
            end
        end
        (columns=names, data=permutedims(reshape(values, length(names), :)))
    end
end

function component_index(component)
    index = findfirst(==(Symbol(component)), COMPONENTS)
    isnothing(index) && error("Unknown component $component; choose one of $(join(COMPONENTS, ", ")).")
    index
end

"""Return a time × radius × depth view of a field. Editing it edits result.data."""
field(result, component) = selectdim(result.data, 4, component_index(component))

"""Get a radius × depth × component snapshot at a physical time in seconds.

An exact stored time returns a view; editing it edits result.data. Interior
times are linearly interpolated into a new array, requiring finite endpoints.
Interpolation cannot cross a jump or slope change in the source history;
use sample_field for explicit coordinate queries and saved-time intervals.
"""
function snapshot(result, time::Real)
    times = result.time
    first(times) <= time <= last(times) || error("Snapshot time must lie within the saved time window.")
    i = searchsortedfirst(times, time)
    if times[i] == time
        plane = view(result.data, i, :, :, :)
    else
        _check_query_time_interval(result.input,times[i-1],times[i])
        lower, upper = view(result.data, i-1, :, :, :), view(result.data, i, :, :, :)
        all(isfinite, lower) && all(isfinite, upper) ||
            error("Cannot interpolate a snapshot across nonfinite samples. Use a stored time or a physically finite source rise time.")
        fraction = (time-times[i-1])/(times[i]-times[i-1])
        plane = (1-fraction).*lower .+ fraction.*upper
    end
    (time=Float64(time), r=result.input.r, z=result.input.zr,
        components=COMPONENTS, units=RESULT_UNITS, data=plane)
end
