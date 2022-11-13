using PressureTraceXT
using Test

@testset "PressureTraceXT.jl" begin
    @testset "Load traces" begin
        ptrace_20210708run3 = PressureTraceXT.load_wistl_trace(joinpath(@__DIR__, "traces/ptrace_20210708run3.lvm"))
    end
end
