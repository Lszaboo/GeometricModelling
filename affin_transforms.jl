using Juliagebra
using LinearAlgebra

App()

ParametricSurface(callback::Function,width,height,uStart,uEnd,vStart,vEnd,color,dependents::Juliagebra.DependentsT) =
Juliagebra._ParametricSurface(_call=callback,_width=width,_height=height,_uStart=uStart,_uEnd=uEnd,_vStart=vStart,_vEnd=vEnd,_color=color,_deps=dependents)

# ? Helper Triangle function
function Triangle(a,b,c,color)
    return ParametricSurface(3,3,0.0,1.0,0.0,1.0,color,[a,b,c]) do u,v,a,b,c
        
        if (u>=0.5 && v>=0.5)
            u = 0.5
            v = 0.5
        end

        return (1-u-v) .* a[:xyz] .+ u .* b[:xyz] .+ v .* c[:xyz]
    end
end

# ? Helper Corner function
function Corner(o,x,y)
    pntX = x._x + y._x - o._x
    pntY = x._y + y._y - o._y
    pntZ = x._z + y._z - o._z
    
    return Point(pntX,pntY,pntZ) 
end

# ? Helper Quad function
function Quad(p1,p2,p3,p4,color)
    Triangle(p1,p2,p3,color)
    Triangle(p4,p3,p2,color)
end

# ? Helper Cuboid function
function assembleCuboidFaces(p1,p2,p3,p4,p5,p6,p7,p8,color)
    Quad(p1,p3,p2,p5,color)
    Quad(p1,p2,p4,p6,color)
    Quad(p1,p4,p3,p7,color)
    Quad(p4,p6,p7,p8,color)
    Quad(p3,p7,p5,p8,color)
    Quad(p5,p8,p2,p6,color)
end

# ? Helper Cuboid function
function Cuboid(x,y,z,color)
    bottomLeftCorner = (x,y,z) .- 0.5
    
    p1 = Point(bottomLeftCorner...)
    p2 = Point((bottomLeftCorner .+ (1,0,0))...)
    p3 = Point((bottomLeftCorner .+ (0,1,0))...)
    p4 = Point((bottomLeftCorner .+ (0,0,1))...)

    p5 = Corner(p1,p2,p3)
    p6 = Corner(p1,p2,p4)
    p7 = Corner(p1,p3,p4)
    p8 = Corner(p4,p6,p7)

    assembleCuboidFaces(p1,p2,p3,p4,p5,p6,p7,p8,color)

    return [p1,p2,p3,p4,p5,p6,p7,p8]
end

function Sphere(r,centerX,centerY,centerZ,color,toggle)
    return ParametricSurface(5,5,0.0,pi,0.0,2.0*pi,color,[toggle]) do alfa,beta,toggle
        if (toggle[:state])
            return nothing
        end
        
        x = centerX + r * sin(alfa) * cos(beta)
        y = centerY + r * sin(alfa) * sin(beta)
        z = centerZ + r * cos(alfa)

        return (x,y,z)
    end
end

function PlanesOfXYZ()
    # ? Toggle to show xyPlane
    xyToggle = Toggle()
    xyPlane = ParametricSurface(2,2,-10.0,10.0,-10.0,10.0,(1.0,1.0,1.0) .* 3.0,[xyToggle]) do u,v,xyToggle
        if(!xyToggle[:state])
            return nothing
        end

        return (u,v,0.0)
    end

    # ? Toggle to show xzPlane
    xzToggle = Toggle()
    xzPlane = ParametricSurface(2,2,-10.0,10.0,-10.0,10.0,(1.0,1.0,1.0) .* 3.0,[xzToggle]) do u,v,xzToggle
        if(!xzToggle[:state])
            return nothing
        end

        return (u,0.0,v)
    end

    # ? Toggle to show yzPlane
    yzToggle = Toggle()
    yzPlane = ParametricSurface(2,2,-10.0,10.0,-10.0,10.0,(1.0,1.0,1.0) .* 3.0,[yzToggle]) do u,v,yzToggle
        if(!yzToggle[:state])
            return nothing
        end

        return (0.0,u,v)
    end
end

function LinesOfXYZ()
    xToggle = Toggle()
    for x in -10:10
        Sphere(0.125,x,0.0,0.0,(8.0,0.0,0.0),xToggle)
    end
    
    yToggle = Toggle()
    for y in -10:10
        Sphere(0.125,0.0,y,0.0,(0.0,8.0,0.0),yToggle)
    end

    zToggle = Toggle()
    for z in -10:10
        Sphere(0.125,0.0,0.0,z,(0.0,0.0,8.0),zToggle)
    end
    
    

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
        0.0 sin(alfa)  cos(alfa) 0.0;
        0.0 0.0 0.0 1.0
    ]
end

function Rotate010(alfa)
    return [
        cos(alfa) 0.0 sin(alfa) 0.0;
        0.0 1.0 0.0 0.0;
        -sin(alfa) 0.0 cos(alfa) 0.0;
        0.0 0.0 0.0 1.0
    ]
end

function Rotate001(alfa)
    return [
        cos(alfa) -sin(alfa) 0.0 0.0;
        sin(alfa) cos(alfa) 0.0 0.0;
        0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 1.0
    ]
end

function Scale(x,y,z)
    return [
        x 0 0 0;
        0 y 0 0;
        0 0 z 0;
        0 0 0 1.0
    ]    
end

MirrorXY() = Scale(1.0,1.0,-1.0)
MirrorXZ() = Scale(1.0,-1.0,1.0)
MirrorYZ() = Scale(-1.0,1.0,1.0)

LinesOfXYZ()
PlanesOfXYZ()

# ? Toggle for visualizing transformed object
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





cuboidCorners = Cuboid(0,0,0,(0.0,1.0,0.0))
transformedCorners = []

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

    push!(transformedCorners,transformedPoint)
end

assembleCuboidFaces(transformedCorners...,(0.0,0.0,1.0))

play!()