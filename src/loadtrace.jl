# Load WiSTL trace files of the format
#=
LabVIEW Measurement	
Writer_Version	2
Reader_Version	2
Separator	Tab
Decimal_Separator	.
Multi_Headings	Yes
X_Columns	No
Time_Pref	Absolute
Operator	aames
Date	2021/07/08
Time	11:59:39.8490078999995408735
***End_of_Header***	
	
Channels	12												
Samples	40000	40000	40000	40000	40000	40000	40000	40000	40000	40000	40000	40000	
Date	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	2021/07/08	
Time	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	11:59:39.8490078999995408735	
Y_Unit_Label	Volts	Volts	Volts	Volts	Volts	Volts	Volts	Volts	Volts	Volts	Volts	Volts	
X_Dimension	Time	Time	Time	Time	Time	Time	Time	Time	Time	Time	Time	Time	
X0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	0.0000000000000000E+0	
Delta_X	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	1.000000E-6	
***End_of_Header***													
X_Value	PT 1	PT 2	PT 3	PT 4	PT 5	PT 6	PT 7	PT 8	PT 9	PT 10	PT 11	PT 12	Comment
	0.299603	0.238652	0.534667	-0.122070	0.301557	0.534771	0.237630	0.610352	0.124918	-0.061035	0.000000	0.244141
	0.299603	0.536966	0.950518	0.122070	0.603114	0.831865	0.356444	0.610352	0.374755	-0.061035	0.302394	0.732422
[...]
=#
function load_wistl_trace(lvm_path)
    open(lvm_path) do f
        meta = Dict{String,String}(Pair(split(readline(f), '\t', limit = 3)[1:2]...) for i = 1:22)
        header = split(readline(f), '\t')[1:(parse(Int, meta["Channels"]) + 1)] .|> String

        seekstart(f)
        ptrace = CSV.File(f; header, delim='\t', skipto=24, normalizenames=true, drop=[1,]) |> DataFrame

        ptrace, meta
    end
end

"""
    reorder_traces!(ptrace, channelmap_path)

    ptrace::DataFrame contains the traces to be reordered/renamed
    channelmap_path::String path to the channelmap csv file with columns "truename", "givenname", and "keep".

The initial column names of `ptrace` should correspond to `givenname` in the channelmap file.
The output order of the columns will be the same as the order of the `truename` column.
"""
function reorder_traces!(ptrace, channelmap_path)
    channelmap = CSV.File(channelmap_path) |> DataFrame
    
    mapping = channelmap[channelmap.keep, :givenname] .=> channelmap[channelmap.keep, :truename]
    select!(ptrace, mapping)
end

"""
    PressureTrace(lvm_path, channelmap_path)

    lvm_path::String path to the WiSTL trace file
    channelmap_path::String path to the channelmap csv file with columns "truename", "givenname", and "keep".

The filepath constructor for `PressureTrace` loads a raw WiSTL trace file, reorders/renames
the columns according to the channelmap file, inserts a column for the time in seconds,
and returns a `PressureTrace`.
"""
struct PressureTrace
    meta::Dict{String, String}
    data::DataFrame
    time::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    function PressureTrace(lvm_path::String, channelmap_path::String)
        ptrace, meta = load_wistl_trace(lvm_path)
        reorder_traces!(ptrace, channelmap_path)        
        time = range(parse(Float64, meta["X0"]), length = size(ptrace, 1), step = parse(Float64, meta["Delta_X"]))
        # insertcols!(ptrace, 1, :time => time)
        # PressureTrace(meta, ptrace)
        new(meta, ptrace, time)
    end
end

function Dates.DateTime(ptrace::PressureTrace)
    date = Date(ptrace.meta["Date"], "yyyy/mm/dd")
    time = Time(first(split(ptrace.meta["Time"], '.')))
    DateTime(date, time)
end

function Base.show(io::IO, ptrace::PressureTrace)
    print(io, "PressureTrace with $(size(ptrace.data, 2)) traces recorded at $(Dates.DateTime(ptrace))")
end