# basis functions ('G functions') in both Euclidean and Minkowski momenta.

using SpecialFunctions

# in Euclidean Momentum
GEpp(c::Float64,p::Float64) = besselk(c+0.5,p*zratio)*besseli(c+0.5,p) - besseli(c+0.5,p*zratio)*besselk(c+0.5,p)
GEpm(c::Float64,p::Float64) = besselk(c+0.5,p*zratio)*besseli(c-0.5,p) + besseli(c+0.5,p*zratio)*besselk(c-0.5,p)
GEmp(c::Float64,p::Float64) = besselk(c-0.5,p*zratio)*besseli(c+0.5,p) + besseli(c-0.5,p*zratio)*besselk(c+0.5,p)
GEmm(c::Float64,p::Float64) = besselk(c-0.5,p*zratio)*besseli(c-0.5,p) - besseli(c-0.5,p*zratio)*besselk(c-0.5,p);

GE11(p::Float64) = besselk(1,p*zratio)*besseli(1,p) - besseli(1,p*zratio)*besselk(1,p)
GE00(p::Float64) = besselk(0,p*zratio)*besseli(0,p) - besseli(0,p*zratio)*besselk(0,p)
GE10(p::Float64) = besselk(1,p*zratio)*besseli(0,p) + besseli(1,p*zratio)*besselk(0,p)
GE01(p::Float64) = besselk(0,p*zratio)*besseli(1,p) + besseli(0,p*zratio)*besselk(1,p)



# in Lorentz Momentum
Gmpz(c::Float64,p::Float64,z::Float64) = pi/2*(besselj(c-0.5,p*z)*bessely(c+0.5,p) - bessely(c-0.5,p*z)*besselj(c+0.5,p))
Gpmz(c::Float64,p::Float64,z::Float64) = pi/2*(besselj(c+0.5,p*z)*bessely(c-0.5,p) - bessely(c+0.5,p*z)*besselj(c-0.5,p))
Gmmz(c::Float64,p::Float64,z::Float64) = pi/2*(besselj(c-0.5,p*z)*bessely(c-0.5,p) - bessely(c-0.5,p*z)*besselj(c-0.5,p))
Gppz(c::Float64,p::Float64,z::Float64) = pi/2*(besselj(c+0.5,p*z)*bessely(c+0.5,p) - bessely(c+0.5,p*z)*besselj(c+0.5,p))

Gmp(c::Float64, p::Float64) = Gmpz(c::Float64, p::Float64, zratio) 
Gpm(c::Float64, p::Float64) = Gpmz(c::Float64, p::Float64, zratio)
Gmm(c::Float64, p::Float64) = Gmmz(c::Float64, p::Float64, zratio)
Gpp(c::Float64, p::Float64) = Gppz(c::Float64, p::Float64, zratio)


# in Lorentz Momentum
G01z(p::Float64,z::Float64) = pi/2*(besselj0(p*z)*bessely1(p) - bessely0(p*z)*besselj1(p))
G10z(p::Float64,z::Float64) = pi/2*(besselj1(p*z)*bessely0(p) - bessely1(p*z)*besselj0(p))
G00z(p::Float64,z::Float64) = pi/2*(besselj0(p*z)*bessely0(p) - bessely0(p*z)*besselj0(p))
G11z(p::Float64,z::Float64) = pi/2*(besselj1(p*z)*bessely1(p) - bessely1(p*z)*besselj1(p))

G01(p::Float64) = G01z(p, zratio) 
G10(p::Float64) = G10z(p, zratio)
G00(p::Float64) = G00z(p, zratio)
G11(p::Float64) = G11z(p, zratio)