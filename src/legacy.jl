
using Dates
using DataFrames
using Unitful
using ImageFiltering: mapwindow
using DSP
using PyThermo
using PyThermo: Chemical
using PyThermo.ShockTube: shockcalc
using Statistics: median, mean
using CSV
using LsqFit
using GLMakie, Colors
# using ColorSchemes # I don't think this is used?
using Printf
using Interpolations

function filtertrace!(ptrace::PressureTrace, filter = (median_size = 15, lowpass_mul = 100, filt = Butterworth(2)))
    freq = 1 / step(ptrace.time)
    lowpass = digitalfilter(Lowpass(freq / filter.lowpass_mul, fs = freq), filter.filt)
    mapcols!(ptrace.data) do s
        s̄ = mapwindow(median, s, filter.median_size)
        ŝ = filtfilt(lowpass, s)
        [median((s[i], s̄[i], ŝ[i])) for i in eachindex(s)]
    end
end

# Shock-fitting
shockpressure(t, p0, ps, ts, R) = (ps - p0) * ((tanh(R * (t - ts)) + 1) / 2) + p0
@. shockpressure(t, params) = shockpressure(t, params...)

function shockfit(t, p, p_trigger, windowsize, i0 = 1)
    # Use data points in neighborhood of trigger to curve-fit
    i_s = findnext(>(p_trigger), p, i0)
    isnothing(i_s) && return nothing
    δ = windowsize ÷ 2
    pdata = p[(i_s-δ):(i_s+δ)]
    tdata = t[(i_s-δ):(i_s+δ)]
    p0 = Float64[pdata[1], pdata[end], t[i_s], 1/mean(diff(tdata))]
    fit = curve_fit(shockpressure, tdata, pdata, p0)
end

function shockdetect(t, p, i_trigger, p_trigger, windowsize = 30)
    ## Assumes trigger PT is neither the first nor last PT
    ## Assumes the last PT is at the endwall

    trigfits = Dict(i_trigger => shockfit(t, p[i_trigger], p_trigger, windowsize))
    trigfits[i_trigger].converged || error("Unable to detect shock")
    # Refine trigger threshold to (p0 + ps)/2
    p_trigger = (trigfits[i_trigger].param[1] + trigfits[i_trigger].param[2]) / 2

    pt_idxs = sort(collect(keys(p)))
    for i in pt_idxs
        trigfits[i] = shockfit(t, p[i], p_trigger, windowsize)
    end

    p_shock = median(fit.param[2] for (idx, fit) in trigfits)
    p_reshock = trigfits[pt_idxs[end]].param[2]
    reshock_trigger = (p_shock + p_reshock) / 2

    reshockfits = empty(trigfits)
    i0 = findfirst(>(reshock_trigger), p[pt_idxs[end]])
    for i in reverse(pt_idxs)
        reshockfit = shockfit(t, p[i], reshock_trigger, windowsize, i0)
        isnothing(reshockfit) && break
        reshockfits[i] = reshockfit
    end

    return trigfits, reshockfits
end


struct PressureTrace{T<:AbstractVector}
    t::T #StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}}
    timestamp::DateTime
    data::DataFrame
    # function PressureTrace(t, timestamp, data)
    #     issorted(t) || ArgumentError("Pressure trace timestamps are")
end

# PressureTrace

function PressureTrace(filepath, filter = (median_size = 15, lowpass_mul = 100, filt = Butterworth(2)))
    open(filepath) do f
        header = Dict(Pair(split(readline(f), '\t', limit = 3)[1:2]...) for i = 1:23)

        t = range(parse(Float64, header["X0"]),
            length = parse(Int, header["Samples"]),
            step = parse(Float64, header["Delta_X"]))

        date = Date(header["Date"], "yyyy/mm/dd")
        time = Time(first(split(header["Time"], '.')))
        timestamp = DateTime(date, time)

        seekstart(f)
        data = CSV.File(f, skipto = 24, drop = [1], delim = '\t',
            header = "PT" .* string.(0:12)) |> DataFrame


        if !isnothing(filter)
            freq = 1 / step(t)
            lowpass = digitalfilter(Lowpass(freq / filter.lowpass_mul, fs = freq), filter.filt)
            mapcols!(data) do s
                s̄ = mapwindow(median, s, filter.median_size)
                ŝ = filtfilt(lowpass, s)
                [median((s[i], s̄[i], ŝ[i])) for i in eachindex(s)]
            end
        end
        PressureTrace(t, timestamp, data)
    end
end

"""
    xt_data(ptrace, ptlocs, driver, driven, trigger, ref_PT)
Calculate the time of shock arrival at each pressure transducer.
Uses the driven gas properties and the Δx/Δt between the trigger
and reference pressure transducers to calculate an approximate Mach
number, which is used to determine driver conditions and shocked 
gas pressure for determination of shock arrival time at each transducer.
    Arguments:
    - `ptrace::PressureTrace`
    - `ptlocs::DataFrame`
    - `driver::Fluid`, e.g. `Species("Helium")`
    - `driven::Fluid`, e.g. `Mixture(["Helium", "Acetone"], zs=[0.95, 0.05])`
    - `trigger::Pair{Symbol, Unitful.Quantity}`, e.g. `:PT3 => 10u"psi"`
    - `ref_PT::Symbol`, e.g. `:PT6`
# Examples
```julia-repl
julia> ptlocs = CSV.File("PT_locations.csv") |> DataFrame
12×3 DataFrame
 Row │ name    x_m       σ_m
     │ String  Float64   Float64
─────┼─────────────────────────────-
   1 │ PT1     0.506412  0.00254
   2 │ PT2     1.36207   0.0035921
   3 │ PT3     1.8415    0.0035921
   4 │ PT4     2.05422   0.0035921
   5 │ PT5     2.33997   0.00299529
   6 │ PT6     3.84492   0.00392726
   7 │ PT7     4.49897   0.00392726
   8 │ PT8     5.10222   0.00392726
   9 │ PT9     5.42607   0.00392726
  10 │ PT10    5.98329   0.00392726
  11 │ PT11    6.50716   0.00392726
  12 │ PT12    6.75481   0.00392726
julia> xt_data(PressureTrace("run3/ptrace.lvm"), ptlocs, 
        Species("N2"), Species("Ar"), :PT3 => 10u"psi", :PT6)
((driver = Species(N2, 298.1 K, 4.275e+06 Pa), 
  driven = Species(Ar, 298.1 K, 1.013e+05 Pa), 
 shocked = Species(Ar, 712.5 K, 6.039e+05 Pa), 
      u2 = 429.3754150350808 m s^-1), 
      12×2 DataFrame
 Row │ x         t_shock
     │ Float64   Float64
─────┼────────────────────
   1 │ 0.506412  0.008039
   2 │ 1.36207   0.009323
   3 │ 1.8415    0.010003
   4 │ 2.05422   0.010307
   5 │ 2.33997   0.010695
   6 │ 3.84492   0.012797
   7 │ 4.49897   0.013711
   8 │ 5.10222   0.014637
   9 │ 5.42607   0.015019
  10 │ 5.98329   0.015796
  11 │ 6.50716   0.016535
  12 │ 6.75481   0.016881)
```
"""
function xt_data(ptrace::PressureTrace, ptlocs::DataFrame, driver::Chemical, driven::Chemical, trigger::Pair, ref_PT::Symbol)
    trigger_PT, trigger_thresh = trigger
    trigger_thresh_psi = ustrip(trigger_thresh |> u"psi")

    # Get PT indices (PT1 => 1, PT9 => 9, etc) within ptlocs
    trigger_idx = parse(Int, string(trigger_PT)[3:end])
    ref_idx = parse(Int, string(ref_PT)[3:end])

    # find trigger time index
    trigger_sensor = >(trigger_thresh_psi)
    trigger_ptrace_idx = findfirst(trigger_sensor, ptrace.data[!, trigger_PT])

    ref_ptrace_idx = findfirst(trigger_sensor, ptrace.data[!, ref_PT])
    Δt = ptrace.t[ref_ptrace_idx] - ptrace.t[trigger_ptrace_idx]
    Δx = ptlocs[ref_idx, :x_m] - ptlocs[trigger_idx, :x_m]

    # calculate Mach number between trigger and reference PTs
    W₀ = (Δx / Δt) * u"m/s"
    M₀ = W₀ / soundspeed(driven)

    # determine gas states
    states = shockcalc(driver, driven, M₀)
    shocked_psi = PyThermo.pressure(states.shocked) |> u"psi" |> ustrip
    reflected_psi = PyThermo.pressure(states.reflected) |> u"psi" |> ustrip

    # t_shock =   [isnothing(i) ? NaN : ptrace.t[i] for i in findfirst.(>(ustrip(trigger_thresh_psi)), eachcol(ptrace.data))]
    # t_reflect = [isnothing(i) ? NaN : ptrace.t[i] for i in findfirst.(>((shocked_psi + reflected_psi/2)), eachcol(ptrace.data))]
    shockfits, reshockfits = shockdetect(ptrace.t, Dict(1:12 .=> eachcol(ptrace.data)), trigger_idx, trigger_thresh_psi, 30)
    t_shock = [isnothing(fit) ? NaN : fit.param[3] for fit in get.(Ref(shockfits), axes(ptrace.data, 2), nothing)]
    t_reflect = [isnothing(fit) ? NaN : fit.param[3] for fit in get.(Ref(reshockfits), axes(ptrace.data, 2), nothing)]

    return (; ptrace, states, shock = (; W₀, M₀), xt = DataFrame(:x => ptlocs[!, :x_m],
        :t_shock => t_shock,
        :t_reflect => t_reflect))
end

"""
An alternate `ptrace` syntax allows for Base Julia-only arguments, namely filepaths, strings, and symbols.
# Examples
```julia-repl
julia> xt_data("run1/ptrace.lvm", "PTlocations.csv", "N2", "Ar", :PT3 => 10, :PT6)
((driver = Species(N2, 298.1 K, 1.633e+06 Pa), 
  driven = Species(Ar, 298.1 K, 1.013e+05 Pa), 
 shocked = Species(Ar, 563.7 K, 4.085e+05 Pa), 
      u2 = 316.0718112278267 m s^-1), 
  12×2 DataFrame
   Row │ x         t_shock
       │ Float64   Float64
  ─────┼────────────────────
     1 │ 0.506412  0.007734
     2 │ 1.36207   0.009218
     3 │ 1.8415    0.010002
     4 │ 2.05422   0.010372
     5 │ 2.33997   0.010854
     6 │ 3.84492   0.013367
     7 │ 4.49897   0.014464
     8 │ 5.10222   0.023391
     9 │ 5.42607   0.016024
    10 │ 5.98329   0.016962
    11 │ 6.50716   0.017844
    12 │ 6.75481   0.018261)
```
"""
function xt_data(ptrace_path, ptloc_path, drivergas, drivengas, trigger, ref_PT;
    ptrace_filter = (median_size = 15, lowpass_mul = 100, filt = Butterworth(2)))

    ptrace = PressureTrace(ptrace_path, ptrace_filter)
    ptlocs = CSV.File(ptloc_path) |> DataFrame
    driver = Species(drivergas)
    driven = Species(drivengas)
    trigger_psi = first(trigger) => last(trigger) * u"psi"
    xt_data(ptrace, ptlocs, driver, driven, trigger_psi, ref_PT)
end

function ptrace_plot(x, t_shock, t_reflect, trace_time, trace_pressure; title = "")
    # x, t_shock, t_reflect, trace.p{NTraces}, trace.t, title

    cscheme = Colors.distinguishable_colors(12, [RGB(0, 0, 0), RGB(1, 1, 1)], dropseed = true)
    f = Figure(resolution = (1200, 800))
    f[1:2, 2] = buttongrid = GridLayout(tellheight = false, tellwidth = true)

    activePT = Node(0)
    allPT = @lift($activePT == 0 ? trues(12) : ($activePT .== 1:12))

    buttons = buttongrid[2:13, 1] = [Button(f, label = "PT$i", height = 23, strokecolor = cscheme[i],
        padding = (5, 5, -30, -30),
        strokewidth = @lift(2 + 4 * (i == $activePT))) for i = 1:12]

    shocktog = buttongrid[2:13, 2] = [Toggle(f, height = 23, active = true, toggleduration = 0.03,
        width = 40, padding = (5, 5, -30, -30))
                                      for i = 1:12]
    reshocktog = buttongrid[2:13, 3] = [Toggle(f, height = 23, active = true, toggleduration = 0.03,
        width = 40, padding = (5, 5, -30, -30))
                                        for i = 1:12]
    buttongrid[1, 1] = Label(f, "Pressure\n transducer", rotation = pi / 4, padding = (-20, -20, -10, 20), valign = :bottom)
    buttongrid[1, 2] = Label(f, "Shock", rotation = pi / 4, padding = (-20, -20, -10, 20), valign = :bottom)
    buttongrid[1, 3] = Label(f, "Reshock", rotation = pi / 4, padding = (-20, -20, -10, 20), valign = :bottom)

    for (i, b) in enumerate(buttons)
        on(b.clicks) do c
            activePT[] = (i == activePT[] ? 0 : i)
        end
    end

    disabletrace = [lift((s, rs) -> !(s[]) & !(rs[]), s.active, rs.active) for (s, rs) in zip(shocktog, reshocktog)]
    xtplot = Axis(f[1, 1])
    xtplot.title = title
    xtplot.ylabel = "x [m]"
    # t_shock = Node(t_shock)
    # t_reshock = Node(t_reflect)
    scatterlines!(xtplot, t_shock, x, markercolor = cscheme, markersize = @lift(5 * (1 .+ $allPT)))
    scatterlines!(xtplot, t_reflect, x, markercolor = cscheme, markersize = @lift(5 * (1 .+ $allPT)))
    hidexdecorations!(xtplot, grid = false)

    ptraces = Axis(f[2, 1])
    for i = 1:12
        lines!(ptraces, trace_time, trace_pressure[i],
            color = @lift(RGBA(cscheme[i], $allPT[i])),
            visible = @lift($allPT[i] && !($(disabletrace[i]))))
        vlines!(ptraces, t_shock[i],
            color = :black, visible = @lift($activePT == i), linestyle = :dash)
        vlines!(ptraces, t_reflect[i],
            color = :black, visible = @lift($activePT == i), linestyle = :dashdot)
    end
    ptraces.xlabel = "t [s]"
    ptraces.ylabel = "p [psi]"
    linkxaxes!(xtplot, ptraces)
    f
end

function ptrace_plot(ptrace)
    ptrace_plot(ptrace.xt.x,
        ptrace.xt.t_shock,
        ptrace.xt.t_reflect,
        ptrace.ptrace.t,
        collect(eachcol(ptrace.ptrace.data)),
        title = string(ptrace.ptrace.timestamp, @sprintf("\t M = %0.2f", ptrace.shock.M₀)))
end
