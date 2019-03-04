include("../so5.jl");
include("../../eeff.jl");


# functions from "../ee_to_ff.jl" signatures

# function p2green(p, UW, UVBC, IRBC, UVblkt, qf1, qf2, af1, af2)
#     """
#     p^2*amplitude in RS. 
#     p        = dimensionless center of mass energy, in the unit of zR.
#     UW       = UV rotation * Wilson line 
#     qf1, qf2 = a list of quantum number*metric*(fermion wavefunction)^2 as a function of dimensionless z in the units of zR
#     af1, af2 = fermion boundary kinetic term. 
#     """

# function p2sm(p, Q1, TL1, Q2, TL2)
#     """
#     from SU(2) and eletric charge to p^2*amplitude in Standard Model. 
#     Note that p = dimensionful center of mass energy, in TeV
#     """


# standard model values
function sm(p, initial, final)
    if initial == "eL"
        Q1 = -1; TL1 = -1/2;
    elseif initial == "eR"
        Q1 = -1; TL1 = 0.;
    else
        error("invalid initial state")
    end

    if final == "bL"
        Q2 = -1/3; TL2 = -1/2;
    elseif final == "bR"
        Q2 = -1/3; TL2 = 0.;
    elseif final == "tL"
        Q2 = 2/3; TL2 = 1/2;
    elseif final == "tR"
        Q2 = 2/3; TL2 = 0.;
    elseif final == "Zh"
        Q2 = 0.; TL2 = 1.;
    else 
        error("invalid initial state")
    end
    return p2sm(p, Q1, TL1, Q2, TL2)
end



# helper functions for relations between variables

make_sb(aL, aB) = sqrt(gy^2/g^2*(kpir+aB)/(kpir+aL))

function aL_from_kR_minphi(kR, minphi) 
    LW0 = 0.5*sin(2*minphi)^2*(2/(1-zratio^2))/(mw^2/kR^2)
    LW = LW0*(1+(mw^2/kR^2)/4*(3/2-1/LW0))
    return LW - kpir
end

minphi_from_v2f2(tuning) = asin(sqrt(tuning))/2

function make_sb_from_kR_minphi(kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL
    return make_sb(aL, aB)
end


function make_UW(phi, sb)
    Ub = Matrix{Float64}(I,4,4); UW = Matrix{Float64}(I,4,4);
    cb = sqrt(1-sb^2);
    sw = sin(2*phi)
    cw = cos(2*phi)
    Ub[2,2] = cb; Ub[2,4] = -sb; Ub[4,2] = sb; Ub[4,4] = cb;
    UW[1,1] = (1+cw)/2; UW[1,2] = (1-cw)/2; UW[1,3] = -sw/sqrt(2);
    UW[2,1] = (1-cw)/2; UW[2,2] = (1+cw)/2; UW[2,3] = sw/sqrt(2);
    UW[3,1] = sw/sqrt(2); UW[3,2] = -sw/sqrt(2); UW[3,3] = cw;
    U = *(Ub,UW); 
end


function make_UVblkt(kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL

    return [aL,0.,0.,aB]
end



# generate a function which takes fermion wavefunctions and UV BLKT and calculates p^2*(matrix of green )
function generate_so5p2green(p, kR, minphi)
    pzR = p/kR

    UVBC = [1,0,0,1]
    IRBC = [1,1,0,1]

    UVblkt = make_UVblkt(kR, minphi)
    sb = make_sb_from_kR_minphi(kR, minphi)
    UW = make_UW(minphi, sb)

    so5p2green(qf1, qf2, af1, af2) = 
        p2green(pzR, UW, UVBC, IRBC, UVblkt, qf1, qf2, af1, af2)

    return so5p2green
end


function rs_eeff(p, kR, minphi, f1, f2, arg1, arg2)
    green = generate_so5p2green(p, kR, minphi)
    q1, func1, af1, qf1 = f1(arg1, kR, minphi)
    q2, func2, af2, qf2 = f2(arg2, kR, minphi)
    return green(qf1, qf2, af1, af2)
end


function compare(rs, sm)
    return 100*( (rs./sm).^2 .-1)
end

function reverse_compare(res)
    return res/100 .+1
end

function dev_eeff(initial, final, p, kR, minphi, f1, f2, arg1, arg2)
    return compare(rs_eeff(p, kR, minphi, f1, f2, arg1, arg2), sm(p, initial, final))
end

function dev_eeff(f1, f2, arg1, arg2, p, kR, minphi)
    initial = string(f1)[1:2]
    final = string(f2)[1:2]
    return dev_eeff(initial, final, p, kR, minphi, f1, f2, arg1, arg2)
end


function helicity(c_l, c_r, E, m) # left and right chirality
    k = E > m ? sqrt(E^2-m^2) : 0.
    h_l = (E+k)*c_l + (E-k)*c_r
    h_r = (E-k)*c_l + (E+k)*c_r
    h_m = m*(c_l + c_r) # mixed
    return h_l, h_r, h_m
end

# e.g. f2 = t5
function dev_eett(initial, final, p, kR, minphi, f1, f2, arg1, arg2)

    f2string = string(f2)
    leftstr = f2string[1]*'L'*f2string[2:end]
    rightstr = f2string[1]*'R'*f2string[2:end]
    f2L = eval.(Meta.parse.(leftstr))
    f2R = eval.(Meta.parse.(rightstr))
    
    E = p/2
    m = arg2["m"]

    sm_c_l = sm(p, initial, f2string[1]*'L')
    sm_c_r = sm(p, initial, f2string[1]*'R')
    sm_h_l, sm_h_r, sm_h_m = helicity(sm_c_l, sm_c_r, E, m)

    rs_c_l = rs_eeff(p, kR, minphi, f1, f2L, arg1, arg2)
    rs_c_r = rs_eeff(p, kR, minphi, f1, f2R, arg1, arg2)
    rs_h_l, rs_h_r, rs_h_m = helicity(rs_c_l, rs_c_r, E, m)

    return compare(rs_h_l, sm_h_l), compare(rs_h_r, sm_h_r), compare(rs_h_m, sm_h_m)
end

function dev_eett(f1, f2string, arg1, arg2, p, kR, minphi)
    initial = string(f1)[1:2]
    final = string(f2)
    return dev_eett(initial, final, p, kR, minphi, f1, f2string, arg1, arg2)
end





# normalized left-handed fermion zero-mode wave-function squared and with metric factor
zR_c_factor(c) = abs(c-0.5) < 1e-7 ? kpir : (1-zratio^(1-2*c))/(1-2*c) # zR-normalized
z0_c_factor(c) = abs(c-0.5) < 1e-7 ? kpir : ((1/zratio)^(1-2*c)-1)/(1-2*c) # z0-normalized

f_zero(c,z) = 1/zR_c_factor(c) *z^(-2*c)


# relations among ct, cb, and beta
c_dependency_for_mbmt(ct, cb) = z0_c_factor(-ct) / z0_c_factor(-cb)

function find_beta(ct, cb)
    return atan(sqrt( mb^2/mt^2 / c_dependency_for_mbmt(ct,cb) ))
end

function mbmt(cb, ct, beta)
    return tan(beta)^2 * c_dependency_for_mbmt(ct, cb)
end

function find_cb(ct, beta)
    tempmbmt(cb) = mbmt(cb,ct,beta)
    return nsolve(tempmbmt, mb^2/mt^2 ; lower = ct, upper = 2.)
end



# # top quark UV boundary kinetic term `at`
# function find_at(ct, cb, beta, kR, angle_dependence)
#     term1 = angle_dependence * cos(beta)^2 / (mt/kR*zratio)^2 / z0_c_factor(-ct)
#     ans = term1 - sin(beta)^2*z0_c_factor(cb) - cos(beta)^2*z0_c_factor(ct)
#     # return ans < 0 ? 1e7 : ans
#     return ans
# end 


function update_c(c)
    window = 0.00009 
    shift  = 0.0001
    if abs(c-0.5) < window || abs(c-1.5) < window || abs(c+0.5) < window || abs(c+1.5) < window
        return c + shift
    else
        return c
    end
end

# top quark UV boundary kinetic term `at`



function find_at0(ct, kR, angle_dependence)
    ct = update_c(ct)
    return angle_dependence / (mt/kR*zratio)^2 / z0_c_factor(-ct) - z0_c_factor(ct)
end


function find_at(ct, cb, beta, kR, angle_dependence)
    ct = update_c(ct)
    cb = update_c(cb)
    a0 = angle_dependence / (mt/kR*zratio)^2 / z0_c_factor(-ct) - z0_c_factor(ct)
    p1 = -1/(4*ct+6)*f(2*ct+3)/f(1+2*ct) + 1/(2-4*ct)*f(1-2*ct)/f(1+2*ct)*zratio^(1+2*ct)
    nn = f(1-2*ct)/(1-2*ct) + a0*zratio^(1-2*ct)
    p2 = 1/(1-2*ct)*( f(1+2*ct)/(2+4*ct)*zratio^(1-2*ct) + f(3-2*ct)/(4*ct-6) ) + a0*zratio^(1-2*ct)*( -1/(4*ct+2) + 1/(4*ct-2)*zratio^2 - 1/(4*ct^2-1)*zratio^(1+2*ct) )
    p3 = zratio^(2*cb-2*ct)*f(1-2*cb)/(1-2*cb)+a0*zratio^(1-2*ct)
    dm2 = -mt^2/kR^2 * (p1 + p2/nn) -sin(beta)^2*p3/nn
    afinal = angle_dependence / (mt/kR*zratio)^2 * (1-dm2) / z0_c_factor(-ct) - z0_c_factor(ct)
    # return ans < 0 ? 1e7 : ans
    return afinal
end 

find_at5(ct, cb, beta, kR, minphi) = find_at(ct, cb, beta, kR, 0.5*sin(2*minphi)^2)
find_at4(ct, cb, beta, kR, minphi) = find_at(ct, cb, beta, kR, sin(minphi)^2)


# d["cb"] = find_cb.(d["ct"], beta)


function make_coupling(kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL

    sb = make_sb(aL, aB)

    g5 = g*sqrt(kpir+aL)
    gx = g5*sqrt(sb^2/(1-sb^2))
    
    return [g5, g5, g5, gx]
end

function make_coupling_in_Y(kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL

    sb = make_sb(aL, aB)

    g5 = g*sqrt(kpir+aL)
    gx = g5*sqrt(sb^2/(1-sb^2))
    
    return [g5, g5/sqrt(1-sb^2), g5, g5*sb]
end


function integrate_with_blkt(integrand, a)
    return quadgk(integrand, zratio, 1)[1] + a*zratio*integrand(zratio)
end

function normalization(func, a ; norm_scope = 1:length(func(1.)))
    # return 1 / ( sum(quadgk(func, zratio, 1)[1][norm_scope]) + a*zratio*sum(func(zratio)[norm_scope]) )
    return 1 / sum(integrate_with_blkt(func,a)[norm_scope]) 
end


function make_qf(func, q, kR, minphi ; UV = true)
    gauge_coupling = make_coupling(kR, minphi)

    if UV
        invU = make_UW(-minphi,0)
    else
        invU = 1
    end

    if typeof(q) == Array{Float64,1} # single f and q, already normalized
        return z-> func(z)*invU*(gauge_coupling.*q)
    else
        return z-> invU*(gauge_coupling.*sum(q.*func(z)))
    end
end


function pure_zeromode(c, q, a, kR, minphi)
    temp(z) = f_zero(c,z)
    func(z) = temp(z)/(1+a*zratio*temp(zratio))
    qe = make_qf(func, q, kR, minphi)
    return q, func, a, qe
end








function assign_qn(q)
    T3L = q[1]; T3R = q[2]; T35 = q[3]; X = q[4]; Y = T3R + X; Q = T3L + Y;
    return T3L, T3R, T35, X, Y, Q
end

function assign_qnlist(qlist)
    if length(qlist[1])  == 1
        return assign_qn(qlist)
    else
        temp = length(qlist)
        T3L, T3R, T35, X, Y, Q = zeros(temp),zeros(temp),zeros(temp),zeros(temp),zeros(temp),zeros(temp)
        for i=1:temp
            T3L[i], T3R[i], T35[i], X[i], Y[i], Q[i] = assign_qn(qlist[i])
        end
        return T3L, T3R, T35, X, Y, Q
    end
end


# note ee is always in the extreme UV
function approxeeff(q1, q2list, func, p, kR, minphi ; a = 0., dQ = 1., KK = 1.)
    LW = kpir + aL_from_kR_minphi(kR, minphi)
    LB = LW
    cw2 = g^2/(g^2+gy^2)
    sw2 = 1-cw2

    T3L1, T3R1, T351, X1, Y1, Q1 = assign_qn(q1)
    T3L2, T3R2, T352, X2, Y2, Q2 = assign_qnlist(q2list)

    one1 = integrate_with_blkt(func, a)
    z2 = integrate_with_blkt(z->func(z)*z^2, a)
    z2logz = 2*integrate_with_blkt(z->func(z)*(z^2*log(1/z)), a)

    delta_g = mz^2/kR^2/4*(3/4 - cw2/LW - sw2/LB)
    delta_s2 = cw2*mz^2/kR^2/4*(1/LW - 1/LB)
    delta_Z = mz^2/kR^2/4*(z2 + z2logz)

    g_eff = g*(1+delta_g)
    s_eff2 = sw2*(1+delta_s2)

    s = sin(2*minphi)
    c = cos(2*minphi)
    delta_Q = sum(( s^2/2*(-T3L2+T3R2) + s*c/sqrt(2)*T352 ).*z2)

    KK_W = sum((-one1/LW + z2 + z2logz).*T3L2)/4
    KK_B = sum((-one1/LB + z2 + z2logz).*Y2)/4

    photon = g^2*sw2*Q1*sum(Q2.*one1)

    QZ1 = T3L1-s_eff2*Q1
    QZ2 = sum((one1+delta_Z).*(T3L2-s_eff2*Q2)) + dQ*delta_Q
    Zboson = p^2/(p^2-mz^2)*g_eff^2/cw2*QZ1*QZ2

    contact = KK* p^2/kR^2 * g^2 * (T3L1*KK_W + sw2/cw2*Y1*KK_B)

    # print("photon : ", round(photon ; sigdigits=4))
    # println(" Zboson : ", round(p^2/(p^2-mz^2)*g_eff^2/cw2*QZ1*sum((one1+delta_Z).*(T3L2-s_eff2*Q2)) ; sigdigits=4))
    # println(" contact : ", round(contact ; sigdigits=3))

    return photon + Zboson + contact

end



# # note ee is always in the extreme UV
# function approxeeff(q1, q2, func, p, kR, minphi ; a = 0.)
#     aL = aL_from_kR_minphi(kR, minphi)
#     aB = aL
#     sb = make_sb(aL, aB)
#     LW = kpir + aL_from_kR_minphi(kR, minphi)
#     LB = LW
#     AW = 1/4*(1/LW-1)
#     AB = 1/4*(1/LB-1)
#     cw2 = g^2/(g^2+gy^2)
#     sw2 = 1-cw2
#     s = sin(2*minphi)
#     c = cos(2*minphi)
#     delta_mz2 = mz^2/kR^2/4*(3/2 -cw2/LW -sw2/LB )
#     s_eff2 = sw2*(1+(AW-AB)*cw2*mz^2/kR^2)

#     T3L1, T3R1, T351, X1, Y1, Q1 = assign_qn(q1)
#     T3L2, T3R2, T352, X2, Y2, Q2 = assign_qn(q2)

#     z2moment  = integrate_with_blkt(z->func(z)*z^2, a)
#     Fc = z2moment + integrate_with_blkt(z->func(z)*(2*z^2*log(1/z)), a)
#     onemoment = integrate_with_blkt(func, a)

#     photon = g^2*sw2*Q1*Q2*onemoment

#     g_eff = g*(1+mz^2/kR^2/4*(3/4-cw2/LW-sw2/LB))
#     compositeQ = mz^2/kR^2/4*Fc*(T3L2-s_eff2*Q2) + ( s^2/2*(-T3L2+T3R2) + s*c/sqrt(2)*T352 )*z2moment
#     QZ1 = T3L1-s_eff2*Q1
#     QZ2 = (T3L2-s_eff2*Q2)*onemoment + compositeQ
#     Zboson = p^2/(p^2-mz^2)*g_eff^2/cw2*QZ1*QZ2

#     contact = p^2/kR^2/4 * g^2/cw2 * (cw2*(Fc-onemoment/LW)*T3L1*T3L2 + sw2*(Fc-onemoment/LB)*Y1*Y2)

#     return photon + Zboson + contact

# end



function approxeeff_zeromode(q1, q2, c, p, kR, minphi ; a = 0.)
    return approxeeff(q1, q2, z->f_zero(c,z), p, kR, minphi ; a = a)
end

function dev_approx_eeff(initial, final, p, kR, minphi, f1, f2, arg1, arg2 ; dQ = 1., KK = 1.)
    if initial == "eL"
        q1 = [-1/2, 0, 0, -1/2]
    elseif initial == "eR"
        q1 = [0, 0, 0, -1]
    else
        error("invalid initial state")
    end

    q2, func2, af2, qf2 = f2(arg2, kR, minphi)

    return compare(approxeeff(q1, q2, func2, p, kR, minphi ; a = af2, dQ = dQ, KK = KK), sm(p, initial, final))
end

function dev_approx_eeff(f1, f2, arg1, arg2, p, kR, minphi ; dQ = 1., KK = 1.)
    initial = string(f1)[1:2]
    final = string(f2)[1:2]
    return dev_approx_eeff(initial, final, p, kR, minphi, f1, f2, arg1, arg2 ; dQ = dQ, KK = KK)
end


function dev_approx_eett(initial, final, p, kR, minphi, f1, f2, arg1, arg2; dQ = 1., KK = 1.)
    if initial == "eL"
        q1 = [-1/2, 0, 0, -1/2]
    elseif initial == "eR"
        q1 = [0, 0, 0, -1]
    else
        error("invalid initial state")
    end

    f2string = string(f2)
    leftstr = f2string[1]*'L'*f2string[2:end]
    rightstr = f2string[1]*'R'*f2string[2:end]
    f2L = eval.(Meta.parse.(leftstr))
    f2R = eval.(Meta.parse.(rightstr))

    E = p/2
    m = arg2["m"]

    sm_c_l = sm(p, initial, f2string[1]*'L')
    sm_c_r = sm(p, initial, f2string[1]*'R')
    sm_h_l, sm_h_r, sm_h_m = helicity(sm_c_l, sm_c_r, E, m)

    q2L, func2L, af2L, qf2L = f2L(arg2, kR, minphi)
    q2R, func2R, af2R, qf2R = f2R(arg2, kR, minphi)
    rs_c_l = approxeeff(q1, q2L, func2L, p, kR, minphi ; a = af2L, dQ = dQ, KK = KK)
    rs_c_r = approxeeff(q1, q2R, func2R, p, kR, minphi ; a = af2R, dQ = dQ, KK = KK)
    rs_h_l, rs_h_r, rs_h_m = helicity(rs_c_l, rs_c_r, E, m)

    return compare(rs_h_l, sm_h_l), compare(rs_h_r, sm_h_r), compare(rs_h_m, sm_h_m)
end

function dev_approx_eett(f1, f2, arg1, arg2, p, kR, minphi; dQ = 1., KK = 1.)
    initial = string(f1)[1:2]
    final = string(f2)
    return dev_approx_eett(initial, final, p, kR, minphi, f1, f2, arg1, arg2; dQ = dQ, KK = KK)
end





function evaluate(dev_func, initial, final, p, kR, minphi, f1, f2, arg1, arg2)
    return dev_func(initial, final, p, kR, minphi, f1, f2, arg1, arg2)
end

function evaluate(dev_func, f1, f2, arg1, arg2, p, kR, minphi)
    initial = string(f1)[1:2]
    final = string(f2)[1:2]
    return dev_func(initial, final, p, kR, minphi, f1, f2, arg1, arg2)
end



















function eL5(args, kR, minphi)
    # fix electron part. assume left-handed e
    ce = 3.
    ae = 0.

    # only UW version
    q = [-0.5, -0.5, 0, 0]
    func(z) = f_zero(ce,z)
    qe = make_qf(func, q, kR, minphi)

    # # Ub*UW version
    # gauge_coupling = make_coupling_in_Y(kR, minphi)
    # sb = make_sb_from_kR_minphi(kR, minphi)
    # invU = inv(make_UW(minphi, sb))
    # qe = [-0.5, -0.5*(1-sb^2), 0, -0.5]
    # # qe = [-0.5, 0, 0, -0.5]
    # qfe(z)  = f_zero(ce, z)*invU*(gauge_coupling.*qe)

    return q, func, ae, qe
end

function eR5(args, kR, minphi)
    # fix electron part. assume left-handed e
    ce = -3.
    ae = 0.

    # only UW version
    q = [0, 0, 0, -1.]
    func(z) = f_zero(-ce,z)
    qe = make_qf(func, q, kR, minphi)

    return q, func, ae, qe
end


function bL5(args, kR, minphi)

    ct = args["ct"]
    cb = args["cb"]
    beta = args["beta"]
    ab = args["ab"]

    cbeta = cos(beta)
    sbeta = sin(beta)

    ctcbratio = f_zero(ct,zratio) / f_zero(cb,zratio)

    func(z) = [cbeta^2 * f_zero(ct,z), # in ct multiplet
             sbeta^2 *ctcbratio* f_zero(cb,z)] # in cb multiplet
    pre = normalization(func, ab)

    q = pre*[  [-1/2, -1/2, 0, 2/3],
            [-1/2, 1/2, 0, -1/3]]

    qb = make_qf(func, q, kR, minphi)

    return q, func, ab, qb
end

function bL5approx(args, kR, minphi)

    ct = args["ct"]
    ab = args["ab"]

    func(z) = f_zero(ct,z) # in ct multiplet
    pre = normalization(func, ab ; norm_scope = 1)

    q = pre*[-1/2, -1/2, 0, 2/3]
    qb = make_qf(func, q, kR, minphi)

    return q, func, ab, qb
end

function bR5(args, kR, minphi)
    cb = args["cb"]
    ab = 0.
    q = [0, 0, 0, -1/3]
    func(z) = f_zero(-cb,z)
    # qb = make_qf(z->f_zero(-cb, z), q, ab, kR, minphi ; UV = false)
    qb = make_qf(func, q, kR, minphi)
    return q, func, ab, qb
end



function t5(args, kR, minphi; t_func = tL_func)
    
    ct = args["ct"]
    cb = args["cb"]
    beta = args["beta"]
    at = args["at"]

    # cbeta = 1.
    # sbeta = 0.
    cbeta = cos(beta)
    sbeta = sin(beta)

    cw = cos(2*minphi)
    sw = sin(2*minphi)
    mtzR = args["m"]/kR

    tpm, tpp, left = t_func(mtzR)
    at = left ? at : 0.
    ctcbratio = ( Gpp(ct, mtzR)*Gpm(ct, mtzR) ) / ( Gpp(cb, mtzR)*Gpm(cb, mtzR) )

    func(z) = [ (-sbeta * ctcbratio * tpm(cb,z))^2,
              (cbeta * (1+cw)/2 * tpm(ct,z))^2,
              (cbeta * (1-cw)/2 * tpm(ct,z))^2,
              (cbeta * (-sw)/sqrt(2) * tpp(ct,z))^2,
              cbeta^2 * (-sw)/sqrt(2) * tpp(ct,z) * ((1+cw)/2 + (1-cw)/2) * tpm(ct,z) ]
    pre = normalization(func, at ; norm_scope = 1:4)

    q = pre*[[0.5, 0.5, 0, -1/3],
            [0.5, -0.5, 0, 2/3],
            [-0.5, 0.5, 0, 2/3],
            [0, 0, 0, 2/3],
            [0, 0, 1, 0]          ]

    qt = make_qf(func, q, kR, minphi ; UV = false)

    UW = make_UW(minphi,0.)
    return [UW*qq for qq in q], func, at, qt

end

function tL_func(mtzR)
    tpm(c,z) = mtzR * sqrt(z) * Gpp(c, mtzR) * Gpmz(c, mtzR, z)
    tpp(c,z) = mtzR * sqrt(z) * Gpm(c, mtzR) * Gppz(c, mtzR, z)
    return tpm, tpp, true
end

function tR_func(mtzR)
    tpm(c,z) = mtzR * sqrt(z) * Gpp(c, mtzR) * Gmmz(c, mtzR, z)
    tpp(c,z) = mtzR * sqrt(z) * Gpm(c, mtzR) * Gmpz(c, mtzR, z)
    return tpm, tpp, false
end

# left-chirality
function tL5(args, kR, minphi)
    return t5(args, kR, minphi ; t_func = tL_func)
end

function tR5(args, kR, minphi)
    return t5(args, kR, minphi ; t_func = tR_func)
end










function Bpp(c)
    c = update_c(c)
    return ( -f(2*c+3)/(4*c+6) + zratio^(1+2*c)*f(1-2*c)/(2-4*c) ) / f(1+2*c)
end

function Bmm(c)
    return Bpp(-c)
end

function Bpm(c,z)
    c = update_c(c)
    return -1/(2+4*c) -z^2/(2-4*c) +z^(1+2*c)/(1-4*c^2)
end

function Bpm(c)
    return Bpm(c,zratio)
end

function Bmp(c,z)
    return Bpm(-c,z)
end

function Bmp(c)
    return Bmp(c,zratio)
end

function La(c,a)
    c = update_c(c)
    return f(1-2*c)/(1-2*c) + a*zratio^(1-2*c)
end

function R(c)
    c = update_c(c)
    return f(1+2*c)/(1+2*c)
end


function Bamm(c, a)
    c = update_c(c)
    return ( f(1-2*c)/(1-2*c)*Bmm(c) + a*zratio^(1-2*c)*Bpm(c) ) / La(c,a)
end


function tL_func_approx(args, kR, angle1, angle2)
    ct = update_c(args["ct"])
    at = args["at"]
    m = args["m"]

    dmt2 =  m^2/kR^2*(-Bpp(ct)-Bamm(ct,at))
    coeff = angle1 + dmt2/2 + m^2/kR^2/2*(-Bpm(ct)-Bamm(ct,at))

    # function func(z)
    #     pre = 1/La(ct,at)*z^(-2*ct)
    #     func1 = coeff*(1+m^2/kR^2*Bpm(ct,z))
    #     func2 = -angle2*( (1-z^(1+2*ct))/(1-zratio^(1+2*ct)) )
    #     return pre*[func1^2, func2^2, func1*func2]
    # end
    function func(z)
        pre = 1/La(ct,at)*z^(-2*ct)
        func11 = (angle1^2 + dmt2 + m^2/kR^2*(-Bpm(ct)-Bamm(ct,at)))*(1+2*m^2/kR^2*Bpm(ct,z))
        func2 = -angle2*( (1-z^(1+2*ct))/(1-zratio^(1+2*ct)) )
        return pre*[func11, func2^2, func2]
    end
    return func, at
end

function tR_func_approx(args, kR, angle1, angle2)
    ct = update_c(args["ct"])
    at = args["at"]
    m = args["m"]

    dmt2 =  m^2/kR^2*(-Bpp(ct)-Bamm(ct,at))
    coeff = 1 + dmt2/2 + m^2/kR^2/2*(-Bpp(ct)+Bpm(ct))

    # function func(z)
    #     pre = 1/R(ct)*z^(2*ct)
    #     func1 = angle2/La(ct,at)*( (1-z^(1-2*ct))/(1-2*ct) )
    #     func2 = coeff*(1+m^2/kR^2*Bmp(ct,z))
    #     return pre*[func1^2, func2^2, func1*func2]
    # end
    function func(z)
        pre = 1/R(ct)*z^(2*ct)
        func1 = angle2/La(ct,at)*( (1-z^(1-2*ct))/(1-2*ct) )
        func22 = (1 + dmt2 + m^2/kR^2*(-Bpp(ct)+Bpm(ct)))*(1+2*m^2/kR^2*Bmp(ct,z))
        return pre*[func1^2, func22, func1]
    end
    return func, 0.
end


function t5approx(args, kR, minphi; t_func = tL_func_approx)

    c = cos(2*minphi)
    s = sin(2*minphi)
    func, at = t_func(args, kR, (1+c)/2, s/sqrt(2))

    q = [[0.5, -0.5, 0, 2/3],
        [0, 0, 0, 2/3],
        [0, 0, 1, 0]          ]
    qt = make_qf(func, q, kR, minphi ; UV = false)

    UW = make_UW(minphi,0.)
    return [UW*qq for qq in q], func, at, qt
end

function tL5approx(args, kR, minphi)
    return t5approx(args, kR, minphi; t_func = tL_func_approx)
end

function tR5approx(args, kR, minphi)
    return t5approx(args, kR, minphi; t_func = tR_func_approx)
end
















# function approxeeff4(q2, func, p, kR, minphi ; a = 0)
#     return approxeeff([-1/2, 0, 0, -1/2], q2, func, p, kR, minphi ; a = a)
# end

function eL4(args, kR, minphi)
    # fix electron part. assume left-handed e
    ce = 3.
    ae = 0.

    q = [-0.5, 0, 0, -0.5]
    func(z) = f_zero(ce,z)
    qe = make_qf(func, q, kR, minphi)

    return q, func, ae, qe
end

function eR4(args, kR, minphi)
    # fix electron part. assume left-handed e
    ce = -3.
    ae = 0.

    q = [0, 0, 0, -1.]
    func(z) = f_zero(-ce,z)
    qe = make_qf(func, q, kR, minphi)

    return q, func, ae, qe
end




function bL4(args, kR, minphi)

    ct = args["ct"]
    cb = args["cb"]
    beta = args["beta"]
    ab = args["ab"]
    
    cbeta = cos(beta)
    sbeta = sin(beta)

    ctcbratio = f_zero(ct,zratio) / f_zero(cb,zratio)

    func(z) = [cbeta^2 * f_zero(ct,z), # in ct multiplet
             sbeta^2 *ctcbratio* f_zero(cb,z)] # in cb multiplet
    pre = normalization(func, ab)

    q = pre*[  [-1/2, 0, 0, 1/6],
           [-1/2, 0, 0, 1/6]]

    qb = make_qf(func, q, kR, minphi)

    return q, func, ab, qb
end

function bL4approx(args, kR, minphi)

    ct = args["ct"]
    ab = args["ab"]
    
    func(z) = f_zero(ct,z) # in ct multiplet
    pre = normalization(func, ab ; norm_scope = 1)

    q = pre*[-1/2, 0, 0, 1/6]
    qb = make_qf(func, q, kR, minphi)

    return q, func, ab, qb
end

function bR4(args, kR, minphi)
    cb = args["cb"]
    ab = 0.
    q = [0, -1/2, 0, 1/6]
    func(z) = f_zero(-cb,z)
    # qb = make_qf(z->f_zero(-cb, z), q, ab, kR, minphi ; UV = false)
    qb = make_qf(func, q, kR, minphi)
    return q, func, ab, qb
end


function t4(args, kR, minphi; t_func = tL_func)
    
    ct = args["ct"]
    cb = args["cb"]
    beta = args["beta"]
    at = args["at"]

    cbeta = cos(beta)
    sbeta = sin(beta)

    cw = cos(minphi)
    sw = sin(minphi)
    mtzR = args["m"]/kR

    tpm, tpp, left = t_func(mtzR)
    at = left ? at : 0.

    ctcbratio = ( Gpp(ct, mtzR)*Gpm(ct, mtzR) ) / ( Gpp(cb, mtzR)*Gpm(cb, mtzR) )

    func(z) = [ (cbeta * cw * tpm(ct,z))^2,
                (cbeta * (-sw) * tpp(ct,z))^2,
                cbeta^2 * (-sw) * cw * tpp(ct,z) * tpm(ct,z),
                (sbeta * cw    * ctcbratio * tpm(cb,z))^2,
                (sbeta * (-sw) * ctcbratio * tpm(cb,z))^2,
                sbeta^2 * (-sw) * cw * (ctcbratio * tpm(cb,z))^2           ]
    pre = normalization(func, at ; norm_scope = [1,2,4,5])

    q = pre*[[1/2, 0, 0, 1/6],
            [0, 1/2, 0, 1/6],
            [0, 0, 1/sqrt(2), 0],
            [1/2, 0, 0, 1/6],
            [0, 1/2, 0, 1/6],
            [0, 0, 1/sqrt(2), 0]]

    qt = make_qf(func, q, kR, minphi ; UV = false)

    UW = make_UW(minphi,0.)
    return [UW*qq for qq in q], func, at, qt

end

function tL4(args, kR, minphi)
    return t4(args, kR, minphi ; t_func = tL_func)
end

function tR4(args, kR, minphi)
    return t4(args, kR, minphi ; t_func = tR_func)
end




function t4approx(args, kR, minphi; t_func = tL_func_approx)

    c = cos(minphi)
    s = sin(minphi)
    func, at = t_func(args, kR, c, s)

    q = [[1/2, 0, 0, 1/6],
        [0, 1/2, 0, 1/6],
        [0, 0, 1/sqrt(2), 0]]
    qt = make_qf(func, q, kR, minphi ; UV = false)

    UW = make_UW(minphi,0.)
    return [UW*qq for qq in q], func, at, qt
end

function tL4approx(args, kR, minphi)
    return t4approx(args, kR, minphi; t_func = tL_func_approx)
end

function tR4approx(args, kR, minphi)
    return t4approx(args, kR, minphi; t_func = tR_func_approx)
end






























function higgs(z)
    return sqrt(2/(1-zratio^2))*z
end

function onshell_Z(kR, minphi)

    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL

    sb = make_sb(aL, aB)

    mzzR = mz/kR

    # in AL3, Z', A35, B order
    wavefunc(z) = z*[ G01z(mzzR*z,zratio/z) + aL*mzzR*zratio*G11z(mzzR*z,zratio/z),
                      G11z(mzzR*z,zratio/z),
                      G11z(mzzR*z,zratio/z),
                      G01z(mzzR*z,zratio/z) + aB*mzzR*zratio*G11z(mzzR*z,zratio/z) ]

    dzwavefunc(z) = z*[ G00z(mzzR*z,zratio/z) + aL*mzzR*zratio*G10z(mzzR*z,zratio/z),
                      G10z(mzzR*z,zratio/z),
                      G10z(mzzR*z,zratio/z),
                      G00z(mzzR*z,zratio/z) + aB*mzzR*zratio*G10z(mzzR*z,zratio/z) ]

    angle_coeff = make_UW(minphi,sb)[1:4,3] # third column vector of Ub*UW
    coeff = angle_coeff./dzwavefunc(1.)
    tempfunc(z) = coeff.*wavefunc(z)

    normal = sum(quadgk(z->1/z*(tempfunc(z)).^2, zratio, 1)[1]) + sum([aL,0,0,aB].*(tempfunc(zratio).^2))

    fullwavefunc(z) = 1/sqrt(normal)*tempfunc(z)
    fulldzwavefunc(z) = 1/sqrt(normal)*coeff.*dzwavefunc(z)

    return fullwavefunc,fulldzwavefunc

end


function approx_onshell_dZ(kR, minphi)
    sw = sin(2*minphi)
    cw = cos(2*minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL
    LW = kpir + aL
    LB = kpir + aB

    sw2 = gy^2/(g^2+gy^2) # weinberg angles
    cw2 = 1-sw2

    sb = make_sb(aL, aB)
    cb = sqrt(1-sb^2)

    mzzR = mz/kR
    delta_mz2 = mz^2/kR^2/4*(3/2 -cw2/LW -sw2/LB )

    ff(z) = z*[-sw/LW*(log(z)+LW), 
                sw*cb, 
                sqrt(2)*cw*(1+mzzR^2/4*(-1/4+1-z^2+2*zratio^2*log(z/zratio))+delta_mz2),
                sw*sb/LB*(log(z)+LB)]
    return ff

end




function Zh(args, kR, minphi)
    vertex_mat = [ 0  0  1  0;
                   0  0 -1  0;
                  -1  1  0  0;
                   0  0  0  0]

    gauge_coupling = make_coupling(kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL
    sb = make_sb(aL, aB)
    invUb = make_UW(0.,-sb)
    invUW = make_UW(-minphi,0.)

    temp, dZfunc = onshell_Z(kR, minphi)

    qZh(z) = 1/z*higgs(z)*gauge_coupling.*(invUW*vertex_mat*invUb*dZfunc(z))

    return 0., 0., 0., qZh
end






function approx_eeZh(p, kR, minphi)
    aL = aL_from_kR_minphi(kR, minphi)
    aB = aL
    sb = make_sb(aL, aB)
    LW = kpir + aL
    LB = kpir + aB
    AW = 1/4*(1/LW-1)
    AB = 1/4*(1/LB-1)
    cw2 = g^2/(g^2+gy^2)
    sw2 = 1-cw2
    s_eff2 = sw2*(1+(AW-AB)*cw2*mz^2/kR^2)
    delta_mz2 = mz^2/kR^2/4*(3/2 -cw2/LW -sw2/LB )
    
    TL3 = -1/2
    Y = -1/2
    Q = -1

    Zboson = p^2/(p^2-mz^2)/(1-delta_mz2) * g^2/cw2 * (1-(cw2*AW+sw2*AB)*mz^2/kR^2) *
                (TL3-s_eff2*Q) * (cos(2*minphi) - 1/16*mz^2/kR^2 + delta_mz2)
    contact = p^2/kR^2 * g^2/cw2 * (-cw2*AW*TL3 + sw2*AB*Y - 1/16*(TL3-sw2*Q))
        
    return Zboson + contact
end





# function eettL4(c1, c2, beta)
#     @printf "processing ct = %.2f \n" c1
    
#     @time begin
    
#     af = fill(zeros(length(d["kR"]), length(d["tuning"])), 6)
    
#     iter = Iterators.product(1:length(d["kR"]),1:length(d["tuning"]))
    
#     af[1] = [find_at4(c1, c2, beta, minphi_from_v2f2(d["tuning"][j]), d["kR"][i]) for (i,j) in iter]
    
#     function make_ff(i,j)
#         kR = d["kR"][i]
#         tuning = d["tuning"][j]
#         minphi = minphi_from_v2f2(tuning)
        
#         cw = cos(minphi)
#         sw = sin(minphi)
#         mtzR = mt/kR

#         pm(c,z) = mtzR * sqrt(z) * Gpp(c, mtzR) * Gpmz(c, mtzR, z)
#         pp(c,z) = mtzR * sqrt(z) * Gpm(c, mtzR) * Gppz(c, mtzR, z)

#         # ct multiplet
#         c1c2ratio = pm(c2,zratio)/pm(c1,zratio)
#         bf11(z) = (cos(beta) * cw * c1c2ratio * pm(c1,z))^2
#         bf12(z) = (cos(beta) * (-sw) * c1c2ratio * pp(c1,z))^2
#         bf13(z) = -sw*cw*(cos(beta) * c1c2ratio)^2 * pm(c1,z) * pp(c1,z)

#         # cb multiplet
#         bf21(z) = (-sin(beta) * cw * pm(c2,z))^2
#         bf22(z) = (-sin(beta) * (-sw) * pm(c2,z))^2
#         bf23(z) = -sw*cw*(-sin(beta))^2 * pm(c2,z) * pm(c2,z)

#         blkt = pm(c2,zratio)^2
        
#         pre = 1/(quadgk(z -> bf11(z)+bf12(z)+bf21(z)+bf22(z), zratio, 1)[1] + af[1][i,j]*zratio*blkt)
#         ff = [bf11, bf12, bf13, bf21, bf22, bf23]
        
#         return pre, ff
#     end
    
#     qf = [  [1/2, 0, 0, 1/6],
#             [0, 1/2, 0, 1/6],
#             [0, 0, 1/sqrt(2), 0],
#             [1/2, 0, 0, 1/6],
#             [0, 1/2, 0, 1/6],
#             [0, 0, 1/sqrt(2), 0]]
    
#     res = fill(zeros(length(d["p"][d["p"].>0.349])), length(d["kR"]), length(d["tuning"]))
#     for i=1:length(d["kR"]), j=1:length(d["tuning"])
#         pre, ff = make_ff(i,j)

#         tuning = d["tuning"][j]
#         minphi = minphi_from_v2f2(tuning)
#         cw = cos(minphi)

#         for k=1:length(ff)
#             res[i,j] += [eeff4(g, minphi_from_v2f2(d["tuning"][j]), ff[k], af[k][i,j] / (cos(beta) * cw)^2, pre*qf[k]) for g in green_mat[i,j][d["p"].>0.349]]
#         end
#     end

#     sm_res = sm_eettL.(d["p"][d["p"].>0.349])
        
#     end
    
#     return [100*( (res[i,j]./sm_res).^2 -1) for (i,j) in iter]
# end



# function eettR4(c1, c2, beta)
#     @printf "processing ct = %.2f \n" c1
    
#     @time begin
    
#     iter = Iterators.product(1:length(d["kR"]),1:length(d["tuning"]))
    
#     function make_ff(i,j)
#         kR = d["kR"][i]
#         tuning = d["tuning"][j]
#         minphi = minphi_from_v2f2(tuning)
        
#         cw = cos(minphi)
#         sw = sin(minphi)
#         mtzR = mt/kR

#         pm(c,z) = mtzR * sqrt(z) * Gpp(c, mtzR) * Gpmz(c, mtzR, z)
#         pp(c,z) = mtzR * sqrt(z) * Gpm(c, mtzR) * Gppz(c, mtzR, z)
#         mm(c,z) = mtzR * sqrt(z) * Gpp(c, mtzR) * Gmmz(c, mtzR, z)
#         mp(c,z) = mtzR * sqrt(z) * Gpm(c, mtzR) * Gmpz(c, mtzR, z)

#         # ct multiplet
#         c1c2ratio = pm(c2,zratio)/pm(c1,zratio)
#         bf11(z) = (cos(beta) * cw * c1c2ratio * mm(c1,z))^2
#         bf12(z) = (cos(beta) * (-sw) * c1c2ratio * mp(c1,z))^2
#         bf13(z) = -sw*cw*(cos(beta) * c1c2ratio)^2 * mm(c1,z) * mp(c1,z)

#         # cb multiplet
#         bf21(z) = (-sin(beta) * cw * mm(c2,z))^2
#         bf22(z) = (-sin(beta) * (-sw) * mm(c2,z))^2
#         bf23(z) = -sw*cw*(-sin(beta))^2 * mm(c2,z) * mm(c2,z)
        
#         pre = 1/(quadgk(z -> bf11(z)+bf12(z)+bf21(z)+bf22(z), zratio, 1)[1])
#         ff = [bf11, bf12, bf13, bf21, bf22, bf23]
        
#         return pre, ff
#     end
    
#     qf = [  [1/2, 0, 0, 1/6],
#             [0, 1/2, 0, 1/6],
#             [0, 0, 1/sqrt(2), 0],
#             [1/2, 0, 0, 1/6],
#             [0, 1/2, 0, 1/6],
#             [0, 0, 1/sqrt(2), 0]]
    
#     res = fill(zeros(length(d["p"][d["p"].>0.349])), length(d["kR"]), length(d["tuning"]))
#     for i=1:length(d["kR"]), j=1:length(d["tuning"])
#         pre, ff = make_ff(i,j)
#         for k=1:length(ff)
#             res[i,j] += [eeff4(g, minphi_from_v2f2(d["tuning"][j]), ff[k], 0, pre*qf[k]) for g in green_mat[i,j][d["p"].>0.349]]
#         end
#     end

#     sm_res = sm_eettR.(d["p"][d["p"].>0.349])
        
#     end
    
#     return [100*( (res[i,j]./sm_res).^2 -1) for (i,j) in iter]
# end