#=

General platform for processes ee -> ff (or more correctly, f1 f1  -> f2 f2 so any fermions but only s-channel)

Green's function is of the form, G(z1)*A*G(z2) - G(z<)*Gbar(z>), where
    1) G : diagonal matrix of G functions which satisfy IR boundary conditions and 
            (+ = 1) bessel function index for z-dependent part
    2) Gbar : diagonal matrix of G functions with (+) bessel function index for z-dependent part
            and oppoite b.c. for IR. extra minus signs for (+) IR b.c.
            i.e. "G11"=>"G10","G10"=>"-G11","G01"=>"G00","G00"=>"-G01"
    3) A : obtained from matrix inversion from b_{UV}*UW*G * A = b_{UV}*UW*Gbar

=#


include("constants.jl");
include("Gfunctions.jl");


"""
    make_string(UVBC, IRBC)

make a matrix of strings G00, G01, G10, G11, of size length(UVBC)*length(IRBC).
This gives a matrix L = b_{UV}*UW*G without the angle factors. 
That is, L_{ij} = U_{ij}*G_{ij} and this function generates the G and Gbar matrix in strings
Gbar is not with the correct sign yet.
"""
function make_string(UVBC, IRBC)


    G_string = ["$(-UV+1)" for UV in UVBC].*reshape(["$(-IR+1)" for IR in IRBC],1,length(IRBC))

    # for Gbar, IR indices are the opposite
    Gbar_string = ["$(-UV+1)" for UV in UVBC].*reshape(["$IR" for IR in IRBC],1,length(IRBC))
    # note Gbar_string is not yet with the correct sign

    return "G".*G_string, "G".*Gbar_string
end



"""
    blkt!(mat, p, UVBC, IRBC, UVblkt)

Apply the UV boundary localized kinetic term `UVblkt` on `mat` in place. 
"""
function blkt!(mat, p, UVBC, IRBC, UVblkt)
    for (i,UV) in enumerate(UVBC)
        if UV == 0 || isequal(0, UVblkt[i]) continue
        else 
            for (j,IR) in enumerate(IRBC)
                mat[i,j] += p*UVblkt[i]*zratio*eval(Meta.parse("G1$(-IR+1)"))(p)
            end
        end
    end
end

function blkt_bar!(mat, p, UVBC, IRBC, UVblkt)
    blkt!(mat, p, UVBC, -IRBC.+1, UVblkt) # for Gbar, IR indices are the opposite
end


# note G => Gbar with correct signs : "G11"=>"G10","G10"=>"-G11","G01"=>"G00","G00"=>"-G01"
function minus_for_bar(IRBC)
    return Diagonal(-2*IRBC.+1)
end


function Amatrix(p, UW, UVBC, IRBC, UVblkt)
    """
    Calculate A matrix from G*A = Gbar
    Note `G` = b_{UV} UW G b_{IR} 
         `Gbar` = b_{UV} UW Gbar b_{IR} 
    """
    @assert length(UVBC) == length(IRBC) == length(UVblkt)
    
    G_string, Gbar_string = make_string(UVBC, IRBC)

    G_expr  = eval.(Meta.parse.(G_string))
    Gbar_expr = eval.(Meta.parse.(Gbar_string))

    G = [temp(p) for temp in G_expr]
    Gbar = [temp(p) for temp in Gbar_expr]

    blkt!(G, p, UVBC, IRBC, UVblkt)
    blkt_bar!(Gbar, p, UVBC, IRBC, UVblkt)

    Gbar = Gbar*minus_for_bar(IRBC)
    
    return *(inv(G.*UW), Gbar.*UW)
end


"generate G(z) ang Gbar(z)"
function make_zstring(IRBC)
    G_string, Gbar_string = make_string([0], IRBC) # always "1" for z1 index
    G_string = reshape(G_string, length(IRBC))
    Gbar_string = reshape(Gbar_string, length(IRBC))
    
    # note Gbar_string is not yet with the correct sign
    return G_string.*"z", Gbar_string.*"z"
end



func_mult(f,g) = z-> [x*y for (x,y) in zip(f(z), g(z))]


function p2green(p, UW, UVBC, IRBC, UVblkt, qf1, qf2, af1, af2)
    """
    p^2*amplitude in RS. 
    p        = dimensionless center of mass energy, in the unit of zR.
    UW       = UV rotation * Wilson line 
    qf1, qf2 = a list of quantum number*metric*(fermion wavefunction)^2 as a function of dimensionless z in the units of zR
    af1, af2 = fermion boundary kinetic term. 
    """
    
    common_factor = p^3 
    # one p from the propagator factor k*p*zR, and p^2 from explicit multiplication of p^2
    
    Amat = Amatrix(p, UW, UVBC, IRBC, UVblkt)
    G_string, Gbar_string = make_zstring(IRBC)
    Gz = eval.(Meta.parse.(G_string))
    Gbarz = eval.(Meta.parse.(Gbar_string))

    G(z) = [z*Gz[i](p,z) for i=1:length(UVBC)]  # don't forget the extra factor of z here!
    Gbar(z) = [z*Gbarz[i](p,z) for i=1:length(UVBC)]
    
    # explicit for-loop (or comprehension) is much faster than broadcasting
    full_z1 = func_mult(qf1, G)
    full_z2 = func_mult(qf2, G)
    
    full1 = quadgk(full_z1, zratio, 1)[1] + af1*zratio*full_z1(zratio)
    full2 = quadgk(full_z2, zratio, 1)[1] + af2*zratio*full_z2(zratio)
    
    full_z1bar = func_mult(qf1, Gbar) # don't forget the extra factor of z here!
    full_z2bar = func_mult(qf2, Gbar)


    # below is correct if af1 == 0, but I am not sure if it is still correct af1 != 0 .....
    after_z2_integration(z1) = full_z1(z1)'*minus_for_bar(IRBC)*(quadgk(full_z2bar, zratio, z1)[1] + af2*zratio*full_z2bar(zratio)) + 
                                full_z1bar(z1)'*minus_for_bar(IRBC)*(quadgk(full_z2, z1, 1)[1])
    
    fullbar = quadgk(after_z2_integration, zratio, 1)[1] + af1*zratio*after_z2_integration(zratio)

    return common_factor*(full1'*Amat*full2 - fullbar)
end



function p2green(p, UW, UVBC, IRBC, UVblkt, q1, q2, f1, f2, af1, af2)
    qf1(z) = f1(z)*q1
    qf2(z) = f2(z)*q2
    return p2green(p, UW, UVBC, IRBC, UVblkt, qf1, qf2, af1, af2)
end
    

function p2sm(p, Q1, TL1, Q2, TL2)
    """
    from SU(2) and eletric charge to p^2*amplitude in Standard Model. 
    Note that p = dimensionful center of mass energy, in TeV
    """
    sw2 = gy^2/(g^2+gy^2) # weinberg angles
    cw2 = 1-sw2
    e2 = g^2*sw2
    return e2*Q1*Q2 + g^2/cw2*(TL1-sw2*Q1)*(TL2-sw2*Q2)*p^2/(p^2-mz^2)
end

