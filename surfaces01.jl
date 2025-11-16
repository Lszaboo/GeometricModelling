using Juliagebra

App()

ParametricSurface(callback::Function,width,height,uStart,uEnd,vStart,vEnd,color,dependents::Juliagebra.DependentsT) =
Juliagebra._ParametricSurface(_call=callback,_width=width,_height=height,_uStart=uStart,_uEnd=uEnd,_vStart=vStart,_vEnd=vEnd,_color=color,_deps=dependents)

# P00 = Point(0, 0, 0)
# P10 = Point(5, 0, 0)
# P01 = Point(0, 5, 0)
# P11 = Point(5, 5, 5)
# 
# bezierSurface = ParametricSurface(50,50,0,1,0,1,[P00,P10,P01,P11]) do u, v, p00, p10, p01, p11
#     return  ((1-u)*(1-v)) .* p00[:xyz] .+ 
#             (   u *(1-v)) .* p10[:xyz] .+
#             ((1-u)*   v)  .* p01[:xyz] .+
#             (   u *   v)  .* p11[:xyz]
# end



#surface = ParametricSurface(50,50,0,1,0,1) do u, v
#    return (u,v,u^2.0 + v^2.0)
#end
#
#curve = ParametricCurve(0.0,1.0,100) do t
#    return b(t)
#end

function constructCurve(f,color=(rand(),rand(),rand()))
    return ParametricCurve(0.0,1.0,100,color,Vector{Juliagebra.PlanDNA}()) do t
        return f(t)
    end
end

function constructSurface(f,color=(rand(),rand(),rand()))
    return ParametricSurface(50,50,0,1,0,1,color,Vector{Juliagebra.PlanDNA}()) do u,v
        return f(u,v)
    end
end

function b(u)
    return (u, 0.0, sin(8*u))
end

function c(u)
    return (u,1.0,u^2)
end

p=1

B(u) = b(u^p)
C(u) = c(u^(1/p))

function x(u,v)
    return (1-v) .* B(u) .+ v .* C(u)
end

# constructCurve(B)
# constructCurve(C)
# constructSurface(x)

function c1(t)
    return (t,0,0)
    #return (t,0,sin(7*t))
end

function c2(t)
    return (0,t,0)
    #return (0,t,sin(8*t))
end

#Point(0,0,2)

function tr(u,v)
    return c1(u) .+ c2(v) .- c1(0)
end

constructCurve(c1)
constructCurve(c2)
# constructSurface(tr)

function d1(t)
    return (t,1.0,0)    
end

function d2(t)
    return (1.0,t,0)    
end

function rc(u,v)
    return (1-v) .* c1(u) .+ v .* c2(u) 
end

function rd(u,v)
    return (1-u) .* d1(v) .+ u .* d2(v) 
end

function rcd(u,v)
    return  ((1-u)*(1-v)) .* p00 .+ 
            (   u *(1-v)) .* p10 .+
            ((1-u)*   v)  .* p01 .+
            (   u *   v)  .* p11
end

function coon(u,v)
    return rc(u,v) .+ rd(u,v) .- rcd(u,v)
    #return rcd(u,v)
end

constructCurve(d1)
constructCurve(d2)
constructSurface(coon)

play!()