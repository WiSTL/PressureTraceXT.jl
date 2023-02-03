function triggerindex(ptrace::PressureTrace, trigger)
    trig_PT, thresh = trigger
    findfirst(>(thresh), ptrace.data[!, trig_PT])
end

"""
    shock_time(ptrace::PressureTrace, trigger)

    ptrace::PressureTrace
    trigger::Pair{Union{String, Symbol}, Float64} the pressure trace column name and threshold
    
    Returns the time of shock detection relative to the trigger time
"""
function shock_time(ptrace::PressureTrace, trigger)
    i0 = triggerindex(ptrace, trigger)
    triggerindices = map(col -> triggerindex(ptrace, col => last(trigger)), names(ptrace.data))
    all(>(0), diff(triggerindices)) || throw(ArgumentError("Trigger indices are not monotonic"))
    ptrace.time[triggerindices] .- ptrace.time[i0]
end

function load_locations(loc_path)
    locs = CSV.File(loc_path) |> DataFrame
    select!(locs, [:name, :x_m, :σ_m] => ((n, x, σ) -> (name=n, x = x .± σ)) => AsTable)
end

function xt(ptrace::PressureTrace, loc_path, trigger)
    xt = load_locations(loc_path)
    filter!(r -> r.name in names(ptrace.data), xt)
    xt.t = shock_time(ptrace, trigger)
    xt
end

"""
stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)

Generates an arbitrarily-spaced finite difference stencil for points x evaluated at x₀
of order m

Ref: https://discourse.julialang.org/t/generating-finite-difference-stencils/85876/5
"""
function stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)
    ℓ = 0:length(x)-1
    m in ℓ || throw(ArgumentError("order $m ∉ $ℓ"))
    A = @. (x' - x₀)^ℓ / factorial(ℓ)
    return A \ (ℓ .== m) # vector of weights w
end

function shockspeed(xt, N=1)
    map(axes(xt, 1)) do i
        window = max(1, i - N):min(size(xt, 1), i + N)
        mapreduce(*, +, stencil(xt.t[window], xt.t[i], 1), xt.x[window])
    end
end
