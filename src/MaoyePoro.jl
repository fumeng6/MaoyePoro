module MaoyePoro

using LinearAlgebra
using Printf
using SpecialFunctions
using FFTW
using TOML
using JLD2

include("types.jl")
include("model.jl")
include("kernel.jl")
include("transform.jl")
include("integration.jl")
include("volume_fluid.jl")
include("adaptive_hankel.jl")
include("response.jl")
include("io.jl")
include("config.jl")
include("archive.jl")
include("queries.jl")
include("tensor_queries.jl")
include("coulomb.jl")
include("cli.jl")

export PoroelasticInput, LayerModel, LayerGrid, read_input, build_model, build_grid,
       solve, run_simulation, write_input, write_outputs, read_timeseries, COMPONENTS,
       save_result, load_result, read_field, field, snapshot, sample_field, strain_tensor,
       read_coulomb_input, run_coulomb, coulomb_from_strain, coulomb_stress, optimal_coulomb_stress

"""Compute displacement, strain, tilt, pore pressure and Darcy flow in SI units.

Array order: time, radius, depth, component. Darcy flow is fluid volume per total
cross-sectional area and time (m/s), rather than pore-space particle velocity.

Source type 1 returns only the pore-pressure-induced correction driven by a
prescribed drained volume change; the original drained elastic field is not included.

`parallel=true` evaluates independent frequencies on the available Julia threads.
It defaults to true when Julia was started with multiple threads.
"""
function solve(input::PoroelasticInput; verbose::Bool=true, max_frequencies=nothing,
               image_factor=DEFAULT_IMAGE_FACTOR, parallel::Bool=Threads.nthreads()>1)
    model = build_model(input)
    grid = build_grid(input, model)
    data, diagnostics = solve_spectra(input, model, grid; verbose, max_frequencies, image_factor, parallel)
    for d in diagnostics
        if !(d.sampling_converged && d.wavenumber_converged && d.frequency_converged)
            @warn "A numerical truncation criterion was not met; inspect result.diagnostics" profile=d.profile sampling=d.sampling_converged wavenumbers=d.wavenumber_converged frequencies=d.frequency_converged
        end
        if d.singular_samples
            @warn "An ideal source-history jump gives a singular instantaneous field; signed Inf samples represent that physical limit" profile=d.profile
        end
    end
    enforce_surface!(data, input, model, grid;effective_radius=first(diagnostics).effective_radius)
    equivalent = nothing
    equivalent_diagnostics = nothing
    if input.compute_equivalent && input.istype != 0
        equivalent, equivalent_diagnostics = equivalent_source(input, model, grid; verbose, max_frequencies, image_factor, parallel)
        direct_ok = !hasproperty(equivalent_diagnostics, :direct_limit_converged) || equivalent_diagnostics.direct_limit_converged
        if !(equivalent_diagnostics.frequency_converged && equivalent_diagnostics.wavenumber_converged && direct_ok)
            @warn "Equivalent-source truncation criterion was not met" equivalent_diagnostics
        end
        if hasproperty(equivalent_diagnostics, :impulsive_steps) && equivalent_diagnostics.impulsive_steps
            @warn "A sudden source change requires an instantaneous fluid-volume pulse; finite rate samples exclude this pulse, and equivalent_diagnostics.impulses records its time and volume."
        end
    end
    if input.istype == 1 && (input.sinj0[1] != 0 || any(input.ts0[i] == input.ts0[i-1] && input.sinj0[i] != input.sinj0[i-1] for i in 2:length(input.ts0)))
        @warn "A sudden prescribed drained volume change has an equivalent instantaneous fluid-volume pulse; finite time samples do not include that pulse volume."
    end
    settings = (max_frequencies=frequency_capacity(input, max_frequencies),
        image_factor=Float64(image_factor), parallel=parallel && Threads.nthreads()>1 &&
            !(input.istype == 0 && all(iszero,model.alf)),
        julia_threads=Threads.nthreads())
    (;input, model, grid, time=collect(range(0.0, input.twindow; length=input.nt)),
       data, diagnostics, equivalent, equivalent_diagnostics, settings)
end

"""Read a TOML input, solve it, and save a complete JLD2 archive by default.

The default output folder is results/<case> beside the input. `output_dir`
overrides it; `format=:text` or `:both` explicitly exports analysis tables.
"""
function run_simulation(path::AbstractString; output_dir=nothing, format=:jld2,
                  verbose::Bool=true, max_frequencies=nothing,image_factor=DEFAULT_IMAGE_FACTOR,
                  parallel::Bool=Threads.nthreads()>1)
    if output_dir === nothing
        output_dir = joinpath(dirname(abspath(path)), "results", splitext(basename(path))[1])
    end
    run_simulation(read_input(path); output_dir, format, verbose, max_frequencies, image_factor, parallel)
end

"""Solve an in-memory input and save one complete JLD2 result by default."""
function run_simulation(input::PoroelasticInput; output_dir::AbstractString="results/simulation", format=:jld2,
                  verbose::Bool=true, max_frequencies=nothing,image_factor=DEFAULT_IMAGE_FACTOR,
                  parallel::Bool=Threads.nthreads()>1)
    Symbol(format) in (:text, :jld2, :both) || error("Output format must be jld2, text, or both.")
    result = solve(input; verbose, max_frequencies, image_factor, parallel)
    paths = write_outputs(result; output_dir, format)
    (;result..., paths)
end

end
