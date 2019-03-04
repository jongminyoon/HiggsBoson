# include ("../packages.jl");

using QuadGK, SpecialFunctions # for so5.jl
using SpecialFunctions # for Gfunctions.jl
using DSP # for smooth in solve.jl
# using DataFrames # for data
using Printf, LinearAlgebra


mutable struct Params
    zR    ::Float64
    ctop  ::Float64
    cvar  ::Float64
    aL    ::Float64
    aB    ::Float64
    af    ::Float64
    minphi::Float64
    mh    ::Float64
    f     ::Float64
    valid ::Bool
end


include("../constants.jl");
include("../Gfunctions.jl");
include("../solve.jl");
# necessary functions : findzR, findf, full_potential, full_ddV, findmh, valid

f(t) = 1-zratio^t;

top_s(t::Float64) = 0.5*sin(2*t)^2
function top_angles(t::Float64) 
    return 0.5*sin(2*t)^2, sin(4*t), 4*cos(4*t), -16*sin(4*t)
end
top_denom(p::Float64, args::Params) = (GEmm(args.ctop, p) + args.af*p*zratio*GEpm(args.ctop, p))*GEpp(args.ctop,p)*(p^2*zratio)


partner_s(t::Float64) = cos(t)^4
function partner_angles(t::Float64) 
    return cos(t)^4, -4*cos(t)^3*sin(t), 3*sin(2*t)^2 -4*cos(t)^4, 6*sin(4*t)+16*cos(t)^3*sin(t)
end
partner_denom(p::Float64, args::Params) = GEmm(args.cvar, p)*GEpp(args.cvar,p)*(p^2*zratio)


boson_s(t::Float64) = 0.5*sin(2*t)^2
function boson_angles(t::Float64) 
    return 0.5*sin(2*t)^2, sin(4*t), 4*cos(4*t), -16*sin(4*t)
end
W_denom(p::Float64, args::Params) = (GE00(p)+args.aL*p*zratio*GE10(p))*GE11(p)*(p^2*zratio)
function Z_denom(p::Float64, args::Params)
    tG00 = GE00(p); tG10 = GE10(p); tG11 = GE11(p);
    sb2 = gy^2/g^2*(kpir+args.aB)/(kpir+args.aL);
    cup2 = (tG00 + args.aB*p*zratio*tG10) / ((tG00 + args.aB*p*zratio*tG10)+sb2*(tG00 + args.aL*p*zratio*tG10))
    temp = cup2*(tG00 + args.aL*p*zratio*tG10)*tG11*(p^2*zratio)
end





function V(p::Float64, t::Float64, args::Params, sine2, denom)
    return log(1 + sine2(t)/denom(p,args))
end

function dV(p::Float64, t::Float64, args::Params, angles, denom)
    s, ds, dds, ddds = angles(t)
    temp = denom(p, args) + s
    return ds/temp
end

function ddV(p::Float64, t::Float64, args::Params, angles, denom)
    s, ds, dds, ddds = angles(t)
    temp = denom(p, args) + s
    return dds/temp - ds^2/temp^2
end

function dddV(p::Float64, t::Float64, args::Params, angles, denom)
    s, ds, dds, ddds = angles(t)
    temp = denom(p, args) + s
    return ddds/temp - 3*dds*ds/temp^2 + 2*ds^3/temp^3
end





function full_V(p::Float64, t::Float64, args::Params)
    return -2*3*V(p, t, args, top_s, top_denom) +
           -2*3*V(p, t, args, partner_s, partner_denom) +
            3*V(p, t, args, boson_s, W_denom) +
            3/2*V(p, t, args, boson_s, Z_denom)
end

function full_dV(p::Float64, t::Float64, args::Params)
    return -2*3*dV(p, t, args, top_angles, top_denom) +
           -2*3*dV(p, t, args, partner_angles, partner_denom) +
            3*dV(p, t, args, boson_angles, W_denom) +
            3/2*dV(p, t, args, boson_angles, Z_denom)
end

function full_ddV(p::Float64, t::Float64, args::Params)
    return -2*3*ddV(p, t, args, top_angles, top_denom) +
           -2*3*ddV(p, t, args, partner_angles, partner_denom) +
            3*ddV(p, t, args, boson_angles, W_denom) +
            3/2*ddV(p, t, args, boson_angles, Z_denom)
end

function full_dddV(p::Float64, t::Float64, args::Params)
    return -2*3*dddV(p, t, args, top_angles, top_denom) +
           -2*3*dddV(p, t, args, partner_angles, partner_denom) +
            3*dddV(p, t, args, boson_angles, W_denom) +
            3/2*dddV(p, t, args, boson_angles, Z_denom)
end





repulsive_denom(p::Float64, args::Params)  = GEmp(args.cvar, p)*GEpm(args.cvar, p)*(p^2*zratio)

function approx_A(arg::Params, denom)
    function temp(x::Float64)
        p = -log(x)
        return p^3/x/denom(p, arg)
    end
    return quadgk(temp,0,1;rtol=1e-10)[1]/(8pi^2)
end














# function findaL(args::Params)
#     temp = mt^2/mw^2*f(1+2*args.ctop)/(1+2*args.ctop)*2/f(2)
#     if abs(args.ctop-0.5) < 1e-4
#         temp *=kpir
#     else
#         temp *=f(1-2*args.ctop)/(1-2*args.ctop)
#     end
#     return temp - kpir
# end

function findaL(args::Params)
    plus_factor  = abs(args.ctop+0.5) < 1e-4 ? kpir : f(1+2*args.ctop)/(1+2*args.ctop)
    minus_factor = abs(args.ctop-0.5) < 1e-4 ? kpir+args.af : f(1-2*args.ctop)/(1-2*args.ctop)+args.af*zratio^(1-2*args.ctop)
    return mt^2/mw^2*plus_factor*minus_factor*2/f(2) - kpir
end

function findzR(args::Params)
    return sqrt(2/f(2)/(kpir+args.aL)/mw^2/2*sin(2*args.minphi)^2)
end

function findf(args::Params)
    return sqrt(vev^2/(1/2*sin(2*args.minphi)^2))
end

function findmh(args::Params)
    if args.f < 1e-3
        return 100.
    else
        temp = integrate(full_ddV, args)
        mh2 = temp(args.minphi)/args.zR^4/(2*args.f^2)
        if mh2 < 0.
            return -1.
        else
            return sqrt(mh2)
        end
    end
end





function printsolve(args::Params)
    # println("cT = ", round(args.cvar,5), ", minphi = ", round(args.minphi,6), ", zR = ", round(args.zR,5), ", mh = ", round(args.mh,5)) 
    temp = args.valid ? "valid" : "invalid"
    Printf.@sprintf("cT = %.5f, minphi = %.5f, zR = %.5f, mh = %.5f, %s \n", args.cvar, args.minphi, args.zR, args.mh, temp)
end;















# TS contour
# PDG values
scenter = 0.02 
sigmas = 0.07
tcenter = 0.06
sigmat = 0.06
corr = 0.92



convmat = [sigmas^2 corr*sigmas*sigmat ; corr*sigmas*sigmat sigmat^2]
major = maximum(sqrt.(eigvals(convmat)))
minor = minimum(sqrt.(eigvals(convmat)))
theta = 0.5*atan(2*corr*sigmas*sigmat/(sigmas^2-sigmat^2))
p68 = sqrt(2.3)
p95 = sqrt(5.991)
t = collect(range(0, step=2*pi/100, length=101))
s68 = major*p68*sin.(t)
t68 = minor*p68*cos.(t)
ss68 = scenter.+([cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[s68' ; t68'])[1,:]
tt68 = tcenter.+([cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[s68' ; t68'])[2,:]

s95 = major*p95*sin.(t)
t95 = minor*p95*cos.(t)
ss95 = scenter.+([cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[s95' ; t95'])[1,:]
tt95 = tcenter.+([cos(theta) -sin(theta) ; sin(theta) cos(theta)]*[s95' ; t95'])[2,:]

mt0 = 0.173
dmt = 0.4
s2 = gy^2/(g^2+gy^2)
c2 = 1-s2
deltaT = 3/8/pi*mt0*dmt/1000/(s2*c2*mz^2)
# Sdmt = zeros(100)
# Tdmt = linspace(-deltaT, +deltaT, 100)
Terr = [-deltaT deltaT]





tth(t::Float64) = cos(2*t)


function findhhh(args::Params)
    if args.f < 1e-3
        return 100.
    else
        temp = integrate(full_dddV, args)
        hhh = temp(args.minphi)/args.zR^4/(sqrt(8)*args.f^3)
        smvalue = 6*2*mh^2/vev
        return hhh/smvalue
    end
end






































