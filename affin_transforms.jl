using Juliagebra
using LinearAlgebra

App()

# ? Helper Triangle function
function Triangle(a,b,c)
    return ParametricSurface(3,3,0.0,1.0,0.0,1.0,[a,b,c]) do u,v,a,b,c
        
        if (u>=0.5 && v>=0.5)
            u = 0.5
            v = 0.5
        end

        return (1-u-v) .* a[:xyz] .+ u .* b[:xyz] .+ v .* c[:xyz]
    end
end

function Quad(o,x,y)
    return ParametricSurface(2,2,0.0,1.0,0.0,1.0,[o,x,y]) do u,v,o,x,y
        
        #if(v>=u)
        #    v=u
        #end
        
        vx = x[:xyz] .- o[:xyz]
        vy = y[:xyz] .- o[:xyz]

        return (o[:xyz] .+ u .* vx .+ v .* vy) 
    end
end

# ? Helper Corner function
function Corner(o,x,y)
    pntX = x._x + y._x - o._x
    pntY = x._y + y._y - o._y
    pntZ = x._z + y._z - o._z

    return Point(pntX,pntY,pntZ,[o,x,y]) do o,x,y
        return  x[:xyz] .+ y[:xyz] .- o[:xyz]
    end
end
# ? Helper Quad function
QuadWithCorner(o,x,y) = return (Quad(o,x,y),Corner(o,x,y))

# ? Helper Corner function
function Corner(o,x,y)
    pntX = x._x + y._x - o._x
    pntY = x._y + y._y - o._y
    pntZ = x._z + y._z - o._z
    
    return Point(pntX,pntY,pntZ) 
end

# ? Helper Quad function

function Quad(p1,p2,p3,p4)
    Triangle(p1,p2,p3)
    Triangle(p4,p3,p2)
end

# ? Helper Cuboid function
function Cuboid(x,y,z)
    bottomLeftCorner = (x,y,z) .- 0.5
    
    p1 = Point(bottomLeftCorner...)
    p2 = Point((bottomLeftCorner .+ (1,0,0))...)
    p3 = Point((bottomLeftCorner .+ (0,1,0))...)
    p4 = Point((bottomLeftCorner .+ (0,0,1))...)

    p5 = Corner(p1,p2,p3)
    p6 = Corner(p1,p2,p4)
    p7 = Corner(p1,p3,p4)
    p8 = Corner(p4,p6,p7)

    Quad(p1,p2,p3,p5)
    Quad(p1,p4,p2,p6)
    Quad(p1,p3,p4,p7)
    Quad(p4,p7,p6,p8)
    Quad(p3,p5,p7,p8)
    Quad(p5,p2,p8,p6)

    return [p1,p2,p3,p4,p5,p6,p7,p8]
end

function Identity()
    return [
        1.0 0 0 0;
        0 1.0 0 0;
        0 0 1.0 0;
        0 0 0 1.0
    ]
end

function Translate(x,y,z)
    return [
        1.0 0 0 x;
        0 1.0 0 y;
        0 0 1.0 z;
        0 0 0 1.0
    ]
end

function Rotate100(alfa)
    return [
        1.0 0.0 0.0 0.0;
        0.0 cos(alfa) -sin(alfa) 0.0;
        0.0 sin(alfa) cos(alfa) 0.0;
        0.0 0.0 0.0 1.0
    ]
end



cuboidCorners = Cuboid(0,0,0)



t1 = Toggle()
txt1 = TextBox("Result = nothing")

gd1 = GenericDependent(Identity(),[txt1]) do txt1
    value = nothing
    
    try
        text = "begin\n" * txt1[:text] * "\nend"
        textSymbols = Meta.parse(text)
    
        UserFunc = @eval function ()
                try 
                $textSymbols
                return Result
            catch err
                println("Error occured evaluating UserFunc:\n$(err)")
            end
            return nothing
        end
        
        value = Base.invokelatest(UserFunc)
    catch err
        println("Error occured parsing UserFunc:\n$(err)\nFor this text:\n$(text)")
    end

    if !(isa(value, AbstractMatrix{Float64}) && size(value) == (4, 4))
        println("Result:\n$(value)\nWasn't a mat4!")
        return Identity()
    end

    return value
end

for corner in cuboidCorners
    transformedPoint = Point(NaN,NaN,NaN,[t1,gd1,corner]) do t1,gd1,corner
        if(!t1[:state])
            return nothing
        end

        xyzw = [corner[:xyz]...,1]

        xyzw = gd1[:val] * xyzw
        w = xyzw[4]
        xyz = tuple(((xyzw / w)[1:3])...)
        
        return xyz 
    end
end


play!()