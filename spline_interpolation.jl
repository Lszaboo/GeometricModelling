using Juliagebra
using JSON
using LinearAlgebra

include("code_box.jl")

App()

struct SplineData
    dataLength::Int
    p::Vector{Tuple{Float64,Float64,Float64}}
    u::Vector{Float64}
    l::Vector{Tuple{Float64,Float64,Float64}}
end

INIT_DATA_LENGTH = 3
INIT_P = [(0.0,0.0,0.0),(1.0,0.0,1.0),(2.0,0.0,-1.0)]
INIT_U = [1.0,2.0,3.0]
INIT_L = [(0.0,0.0,1.0),(0.0,0.0,-1.0),(0.0,0.0,-1.0)]

INIT_SPLINEDATA = SplineData(
    INIT_DATA_LENGTH,
    INIT_P,
    INIT_U,
    INIT_L
)

INIT_SPLINEDATA_STRING = """
{
"U_Modes": ["BOX","UNIFORM","RANDOM","ARCLENGTH"],
"PL_Modes": ["BOX"],
"DataLength": 3,
"p": "BOX",
"u": "BOX",
"l": "BOX"
}
"""

MAX_DATA_LENGTH = 10

NAN_COORDS = (NaN64,NaN64,NaN64)

pData = CodeBox("p",INIT_P)
uData = CodeBox("u",INIT_U)
lData = CodeBox("l",INIT_L)

resetSplineData = Toggle()
txtSplineData = TextBox([resetSplineData]) do resetSplineData
    return INIT_SPLINEDATA_STRING
end


splineData = GenericDependent(INIT_SPLINEDATA,[txtSplineData,pData,uData,lData]) do txtSplineData,pData,uData,lData
    txt = txtSplineData[:text]

    try
        data = JSON.parse(txt)
        
        dataLength = data["DataLength"]

        if (dataLength <=2)
            error("More than 2 points needed!")
        end

        p = []
        pStr = uppercase(data["p"])
        if pStr == "RANDOM"
            for i in 1:dataLength
                xyz = (rand(),rand(),rand()) .* 20.0 .- 10.0 
                push!(p,xyz)
            end
        elseif pStr == "BOX"
            p = pData[:val]
        end
        
        u = []
        uStr = uppercase(data["u"])
        if uStr == "UNIFORM"
            for i in 0:(dataLength-1)
                push!(u,Float64(i))
            end
        elseif uStr == "RANDOM"
            prev = 0.0
            for i in 1:dataLength
                prev = prev + rand()*5.0
                push!(u,prev)
            end
        elseif uStr == "ARCLENGTH"
            
        elseif uStr == "BOX"
            u = uData[:val]
        end
        
        l = []
        lStr = uppercase(data["l"])
        if lStr == "RANDOM"
            for i in 1:dataLength
                xyz = (rand(),rand(),rand()) .* 10.0 .- 5.0 
                push!(l,xyz)
            end
        elseif lStr == "BOX"
            l = lData[:val]
        end

        if !(length(p) == length(u) == length(l) == dataLength)
            error("Length of p,u,v doesn't match DataLength!")
        end

        correctSplineData = SplineData(dataLength,p,u,l)
        println("$(correctSplineData)")
        
        return correctSplineData
    catch e
        println("Error occured:")
        println(string(e))
        println("Reverting to INIT_SPLINEDATA...")
        return INIT_SPLINEDATA
    end

end

points = []
vectorEndPoints = []
uDependents = []

for i in 1:MAX_DATA_LENGTH
    point = Point([splineData]) do splineData
        data = splineData[:val]
        if 0<i && i<=data.dataLength
            return data.p[i]
        end

        return nothing
    end

    vectorEndPoint = Point([splineData]) do splineData
        data = splineData[:val]
        if 0<i && i<=data.dataLength
            return data.p[i] .+ data.l[i]
        end

        return nothing
    end

    uDependent = GenericDependent(0.0,[splineData]) do splineData
        data = splineData[:val]
        if 0<i && i<=data.dataLength
            return data.u[i]
        end

        return NaN64
    end

    Segment(point,vectorEndPoint)

    push!(points,point)
    push!(vectorEndPoints,vectorEndPoint)
    push!(uDependents,uDependent)
end

# ? Hermite-Interpolation

function Bernstein(n,k,t)
    return binomial(n,k) * ((1-t)^(n-k)) * (t^(k))
end

function QuadradicBezier(p0,p1,p2,p3,u,a=0.0,b=1.0)
    u = (u - a) / (b - a)

    val =        p0 .* Bernstein(3,0,u)
    val = val .+ p1 .* Bernstein(3,1,u)
    val = val .+ p2 .* Bernstein(3,2,u)
    val = val .+ p3 .* Bernstein(3,3,u)

    return val
end

function QuadradicBezierCurve(a,b,b0,b1,b2,b3,color)
    return ParametricCurve(0.0,1.0,1000,color,[a,b,b0,b1,b2,b3]) do t,a,b,b0,b1,b2,b3
        u = (b[:val] - a[:val]) * t + a[:val]
        return QuadradicBezier(b0[:val],
                               b1[:val],
                               b2[:val],
                               b3[:val],
                               u,a[:val],b[:val])
    end
end

function CreateHermiteCurve(p0,p1,pl0,pl1,a,b,color)
    b0 = GenericDependent(NAN_COORDS,[p0]) do p0
        return p0[:xyz]
    end

    b1 = GenericDependent(NAN_COORDS,[p0,pl0,a,b]) do p0,pl0,a,b
        mu = b[:val] - a[:val]
        l0 = pl0[:xyz] .- p0[:xyz]
        return l0 .* (mu / 3.0) .+ p0[:xyz]
    end

    b2 = GenericDependent(NAN_COORDS,[p1,pl1,a,b]) do p1,pl1,a,b
        mu = b[:val] - a[:val]
        l1 = pl1[:xyz] .- p1[:xyz]
        return -1 .* l1 .* (mu / 3.0) .+ p1[:xyz]
    end

    b3 = GenericDependent(NAN_COORDS,[p1]) do p1
        return p1[:xyz]
    end
    
    QuadradicBezierCurve(a,b,b0,b1,b2,b3,color)
end

for i in 1:(MAX_DATA_LENGTH-1)
    println("$(i) - $(i+1)")
    p0 = points[i]
    p1 = points[i+1]

    pl0 = vectorEndPoints[i]
    pl1 = vectorEndPoints[i+1]

    a = uDependents[i]
    b = uDependents[i+1]

    CreateHermiteCurve(p0,p1,pl0,pl1,a,b,(0.153,0.773,0.369))    
end

# ? Catmull-Rom-With-Hermite

function FMILLTangent(p0,p2)
    l1 = collect(p2.-p0)
    l1 = normalize(l1)
    return Tuple(l1)
end

function BesselTangentAt0(p0,p1,p2)
    l0 = collect(4.0 .* (p1 .- p0) .- (p2 .- p0))
    l0 = normalize(l0)
    return Tuple(l0)
end

function BesselTangentAt2(p0,p1,p2)
    l2 = collect(4.0 .* (p2 .- p1) .- (p2 .- p0))
    l2 = normalize(l2)
    return Tuple(l2)
end


crLEndPoints = []

p0 = points[1]
p1 = points[2]
p2 = points[3]
crL1 = Point([p0,p1,p2]) do p0,p1,p2
    return p0[:xyz] .+ BesselTangentAt0(p0[:xyz],p1[:xyz],p2[:xyz])
end

push!(crLEndPoints,crL1)

p0 = points[1]
p1 = points[2]
p2 = points[3]
crL2 = Point([p0,p1,p2,splineData]) do p0,p1,p2,splineData
    data = splineData[:val]
    
    if data.dataLength == 2
        return p2[:xyz] .+ BesselTangentAt2(p0[:xyz],p1[:xyz],p2[:xyz])
    end
    
    return p1[:xyz] .+ FMILLTangent(p0[:xyz],p2[:xyz])
end

push!(crLEndPoints,crL2)

for i in 3:(MAX_DATA_LENGTH-1)
    pm1 = points[i-2]
    local p0 = points[i-1]
    local p1 = points[i]
    local p2 = points[i+1]
    crLI = Point([pm1,p0,p1,p2,splineData]) do pm1,p0,p1,p2,splineData
        data = splineData[:val]
        
        if i == data.dataLength
            return p1[:xyz] .+ BesselTangentAt2(pm1[:xyz],p0[:xyz],p1[:xyz])
        elseif i < data.dataLength 
            return p1[:xyz] .+ FMILLTangent(p0[:xyz],p2[:xyz])
        end
        return nothing
    end

    push!(crLEndPoints,crLI)
end

p0 = points[MAX_DATA_LENGTH-2]
p1 = points[MAX_DATA_LENGTH-1]
p2 = points[MAX_DATA_LENGTH]
crLN = Point([p0,p1,p2]) do p0,p1,p2
    return p2[:xyz] .+ BesselTangentAt2(p0[:xyz],p1[:xyz],p2[:xyz])
end

push!(crLEndPoints,crLN)
println(length(crLEndPoints))

for i in 1:(MAX_DATA_LENGTH-1)
    println("$(i) - $(i+1)")
    local p0 = points[i]
    local p1 = points[i+1]

    pl0 = crLEndPoints[i]
    pl1 = crLEndPoints[i+1]

    a = uDependents[i]
    b = uDependents[i+1]

    CreateHermiteCurve(p0,p1,pl0,pl1,a,b,(0.969,0.224,0.439))
end

# ? Catmull-Rom

crB = []

p0 = points[1]
p1 = points[2]
p2 = points[3]
a = uDependents[1]
b = uDependents[2]

crB0 = GenericDependent(NAN_COORDS,[p0,p1,p2,a,b]) do p0,p1,p2,a,b 
    l = BesselTangentAt0(p0[:xyz],p1[:xyz],p2[:xyz])
    mu = b[:val] - a[:val]
    return l .* (mu / 3.0) .+ p0[:xyz]
end

push!(crB,(nothing,crB0))

p0 = points[1]
p1 = points[2]
p2 = points[3]

a = uDependents[1]
b = uDependents[2]
c = uDependents[3]

crB1m1 = GenericDependent(NAN_COORDS,[p0,p1,p2,a,b,c]) do p0,p1,p2,a,b,c
    mum1 = (b[:val] - a[:val]) / (c[:val] - a[:val])
    l = ((p2[:xyz] .- p0[:xyz]) .* mum1) ./ 3.0
    return p1[:xyz] .- l
end

crB1p1 = GenericDependent(NAN_COORDS,[p0,p1,p2,a,b,c]) do p0,p1,p2,a,b,c
    mum1 = (c[:val] - b[:val]) / (c[:val] - a[:val])
    l = ((p2[:xyz] .- p0[:xyz]) .* mum1) ./ 3.0
    return p1[:xyz] .+ l
end

push!(crB,(crB1m1,crB1p1))

for i in 3:(MAX_DATA_LENGTH-1)
    pm1 = points[i-2]
    local p0 = points[i-1]
    local p1 = points[i]
    local p2 = points[i+1]

    local a = uDependents[i-1]
    local b = uDependents[i]
    local c = uDependents[i+1]

    b1m1 = GenericDependent(NAN_COORDS,[pm1,p0,p1,p2,a,b,c,splineData]) do pm1,p0,p1,p2,a,b,c,splineData
        data = splineData[:val]
        
        if i < data.dataLength
            mum1 = (b[:val] - a[:val]) / (c[:val] - a[:val])
            l = ((p2[:xyz] .- p0[:xyz]) .* mum1) ./ 3.0
            return p1[:xyz] .- l
        elseif i == data.dataLength
            l = BesselTangentAt2(pm1[:xyz],p0[:xyz],p1[:xyz])
            mu = b[:val] - a[:val] 
            return -1 .* l .* (mu / 3.0) .+ p1[:xyz]
        end

        return nothing   
    end

    b1p1 = GenericDependent(NAN_COORDS,[p0,p1,p2,a,b,c]) do p0,p1,p2,a,b,c
        mum1 = (c[:val] - b[:val]) / (c[:val] - a[:val])
        l = ((p2[:xyz] .- p0[:xyz]) .* mum1) ./ 3.0
        return p1[:xyz] .+ l
    end

    push!(crB,(b1m1,b1p1))
end

p0 = points[MAX_DATA_LENGTH-2]
p1 = points[MAX_DATA_LENGTH-1]
p2 = points[MAX_DATA_LENGTH]
a = uDependents[MAX_DATA_LENGTH-1]
b = uDependents[MAX_DATA_LENGTH]

crBN = GenericDependent(NAN_COORDS,[p0,p1,p2,a,b]) do p0,p1,p2,a,b
    l = BesselTangentAt2(p0[:xyz],p1[:xyz],p2[:xyz])
    mu = b[:val] - a[:val] 
    return -1 .* l .* (mu / 3.0) .+ p2[:xyz]
end

push!(crB,(crBN,nothing))

for i in 2:MAX_DATA_LENGTH
    local p0 = points[i-1]
    local p1 = points[i]

    local a = uDependents[i-1]
    local b = uDependents[i]

    b0 = GenericDependent(NAN_COORDS,[p0]) do p0
        return p0[:xyz]
    end

    b1 = crB[i-1][2]
    b2 = crB[i][1]

    b3 = GenericDependent(NAN_COORDS,[p1]) do p1
        return p1[:xyz]
    end

    QuadradicBezierCurve(a,b,b0,b1,b2,b3,(0.255,0.706,0.784))
end

play!()