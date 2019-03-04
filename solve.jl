#=

Methods used to find the minimum of Higgs Coleman-Weinberg potential 

=#

using QuadGK



# Euclidean momentum integration
function integrate(integrand, args::Params ; tol = 1e-10)
    """
    integrand is a function which takes three parameters: p, angle t, args

    returns a potential as a function of angle t.
    """
    function (t::Float64)
        function temp(x::Float64)
            p::Float64 = -log(x)
            return p^3/x*integrand(p, t, args)
        end
        return quadgk(temp,0,1;rtol = tol)[1]/(8*π^2)
    end
end



# Find a sub-interval of [x[1], x[5]] which includes the minimum of func(x) in the interval.
# @pre x,y a vector of doubles with length 5
#      y[1,3,5] = func(x[1,3,5]), y[3] < y[5], y[3] < y[1]
function min_five!(x::Vector{Float64}, y::Vector{Float64}, func ; niter::Int64 = 20)
    count::Int64 = 0
    while count != niter
        if x[3] < 1e-6
            x[4] = x[3] + (x[5]-x[3])*0.51
            y[4] = func(x[4])
            if y[3] < y[4]
                x[5] = x[4]; y[5] = y[4];
            else
                x[1] = x[3]; x[3] = x[4];
                y[1] = y[3]; y[3] = y[4];
            end
        else
            x[2] = x[1] + (x[3]-x[1])*0.51
            x[4] = x[3] + (x[5]-x[3])*0.51
            y[2] = func(x[2])
            y[4] = func(x[4])
            if y[2] < y[3] && y[2] < y[4]
                x[5] = x[3]; x[3] = x[2];
                y[5] = y[3]; y[3] = y[2];
            elseif y[3] < y[2] && y[3] < y[4]
                x[1] = x[2]; x[5] = x[4];
                y[1] = y[2]; y[5] = y[4];
            else
                x[1] = x[3]; x[3] = x[4];
                y[1] = y[3]; y[3] = y[4];
            end
        end

        if (x[5]-x[3])/x[3] < 1e-3
            break
        end

        count += 1
    end

end


# initiate x and y for min_five!
function init_five!(x::Vector{Float64}, y::Vector{Float64})

    ind = argmin(y)
    if ind == 1
        x[3] = x[1]; x[5] = x[2];
        y[3] = y[1]; y[5] = y[2];
    elseif ind == 5
        x[3] = x[5]; x[1] = x[4];
        y[3] = y[5]; y[1] = y[4];
    else
        x[5] = x[ind+1]; x[1] = x[ind-1]; x[3] = x[ind];
        y[5] = y[ind+1]; y[1] = y[ind-1]; y[3] = y[ind]; 
    end
end




function valid(args::Params)
    return abs(args.mh - mh) < 0.03/1000.
end

function infeasible(count::Int64, args::Params)
    return count > 9 && (abs((args.mh - mh)/mh) > 0.1 || args.minphi < 0.025)
end

function too_large_minphi(count::Int64, args::Params)
    return count > 9 && args.minphi > 0.3
end



# find the minimum of Higgs potential which is consistent ( valid && not infeasible && not too_large_minphi)
function solve(args::Params)
    upper = 1.
    lower = 0.
    v = integrate(full_V, args)
    t = zeros(Float64,5)
    v_t = zeros(Float64,5)

    small_t = 0.01
    large_t = 0.5
    count = 0

    while true
        args.cvar = (upper + lower)/2.

        t = collect(range(small_t, stop=large_t, length=5))
        v_t = map(v,t)
        init_five!(t, v_t)
        if t[3] > large_t - 0.05
            upper = args.cvar
        else
            min_five!(t, v_t, v)
            args.minphi = t[3]
            args.zR = findzR(args)
            args.f = findf(args)
            args.mh = findmh(args)

            if valid(args)
                args.valid = true
                break
            elseif infeasible(count, args) || too_large_minphi(count, args)
                args.valid = false
                break
            end

            if args.mh < mh
                upper = args.cvar
            else
                lower = args.cvar
            end
        end
        # println(@sprintf("%i th loop, mh= %.5f, minphi = %.5f", count, args.mh, args.minphi))
        count +=1
        if (upper - lower) / (upper + lower) < 1e-6
            args.valid = false
            break
        end
    end
end




# Another solve method utilizing potential derivative and finding dv = 0
function solve1(args::Params)
    v = integrate(full_V, args)
    dv = integrate(full_dV, args)
    t = zeros(Float64,5)
    v_t = zeros(Float64,5)

    small_t = 0.005
    large_t = 0.4

    upper = 1.
    lower = 0.

    args.cvar = upper
    if dv(small_t) > 0
        args.valid = false
        return
    end
    args.cvar = lower
    if dv(large_t) < 0
        args.valid = false
        return
    end

    count = 0
    while true
        args.cvar = (upper + lower)/2.

        if dv(small_t) > 0
            lower = args.cvar
        elseif dv(large_t) < 0
            upper = args.cvar
        else
            t = collect(range(small_t, stop=large_t, length=5))
            v_t = map(v,t)
            init_five!(t, v_t)
            min_five!(t, v_t, v)
            args.minphi = t[3]
            args.zR = findzR(args)
            args.f = findf(args)
            args.mh = findmh(args)

            if valid(args)
                args.valid = true
                break
            elseif infeasible(count, args)
                args.valid = false
                break
            end
            
            if args.mh < mh
                upper = args.cvar
            else
                lower = args.cvar
            end
        end

        if (upper - lower) / (upper + lower) < 1e-6
            args.valid = false
            break
        end
    end
end




# Newtons method for minimization
function newton(f, df, ddf, x0::Float64 ; tol = 1e-10 )
    while true
        dfx = df(x0)
        ddfx = ddf(x0)
        dx = -dfx/ddfx
        lambda2 = dfx^2/ddfx
        if lambda2/2 < tol
            break
        end
        t = backtracking(f, f(x0), dfx, x0, dx)
        x0 += t*dx
    end
    return x0
end

function backtracking(f, fx::Float64, dfx::Float64, x::Float64, dx::Float64 ; alpha = 0.3, beta = 0.8)
    t = 1.
    while f(x+t*dx) > fx + alpha*t*dfx*dx
        t *= beta
    end
    return t
end

function newton_solve(args::Params)
    upper = 1.
    lower = 0.
    args.cvar = (upper + lower)/2.
    v = integrate(full_V, args)
    dv = integrate(full_dV, args)
    ddv = integrate(full_ddV, args)

    small_t = 0.0001
    large_t = 0.4

    for i = 1:25
        print(i, " ")
        if dv(small_t) > 0
            lower = args.cvar
        elseif dv(large_t) < 0
            upper = args.cvar
        else
            args.minphi = newton(v, dv, ddv, small_t)
            args.zR = findzR(args)
            args.f = findf(args)
            args.mh = findmh(args)
            if valid(args)
                break
            end
            if args.mh < mh
                upper = args.cvar
            else
                lower = args.cvar
            end
        end
        args.cvar = (upper + lower)/2.
    end
end










# numerically find solution if the target function is monotonic in its variable
# binary search
function nsolve(func, goal ; lower = 0.5, upper = 1.0, tol = 1e-6)
    c = 0.
    temp = 0.
    upper -= tol
    lower += tol
    count = 0
    
    if isnan(func(upper)) || isnan(func(lower))
        println("NaN found!")
        return nothing
    end
    
    func(upper) > func(lower) ? sign = true : sign = false
    
    for i=1:200
        c = (upper+lower)/2
        temp = func(c)
        (abs(goal) < tol && abs(temp) < tol) && return c # if goal is zero
        abs(temp/goal-1) < tol && return c
        sign ? (temp < goal ? lower = c : upper = c) : (temp > goal ? lower = c : upper = c)
    end

    println("not approaching to the solution")
    return 0.
end







# compute KK masses of various boundary conditions. 
# e.g. pmmass - (+-) bc.


function pmmass(c, zR ; n = 3, af = 0.)
    tempGmp(p::Float64) = Gmp(c,p) + af*p*zratio*Gpp(c,p)
    pp = zeros(n)
    pp[1] = nsolve(tempGmp, 0 ; lower = 0.03, upper = 2.)

    if n > 1
        for i=2:n
            pp[i] = nsolve(tempGmp, 0 ; lower = pp[i-1]+2, upper = pp[i-1]+5)
        end
    end

    return pp./zR
end

function mpmass(c, zR; n = 3)
    tempGpm(p::Float64) = Gpm(c,p)
    pp = zeros(n)
    if c > 0.5
        pp[1] = nsolve(tempGpm, 0 ; lower = 2., upper = 3.5)
    else
        pp[1] = nsolve(tempGpm, 0 ; lower = 1., upper = 2.5)
    end
    if n > 1
        for i=2:n
            pp[i] = nsolve(tempGpm, 0 ; lower = pp[i-1]+2, upper = pp[i-1]+4)
        end
    end

    return pp./zR
end

function ppmass(c, zR; n = 3, af = 0.)
    tempGmm(p::Float64) = Gmm(c,p) + af*p*zratio*Gpm(c,p)
    pp = zeros(n)
    pp[1]= nsolve(tempGmm, 0 ; lower = 1.5, upper = 4.5)
    if n > 1
        for i=2:n
            pp[i] = nsolve(tempGmm, 0 ; lower = pp[i-1]+1.5, upper = pp[i-1]+4)
        end
    end

    return pp./zR
end


function mmmass(c, zR; n = 3)
    tempGpp(p::Float64) = Gpp(c,p)
    pp = zeros(n)
    pp[1] = nsolve(tempGpp, 0 ; lower = 3., upper = 5.)
    if n > 1
        for i=2:n
            pp[i] = nsolve(tempGpp, 0 ; lower = pp[i-1]+2, upper = pp[i-1]+4)
        end
    end

    return pp./zR
end;


function ppbosonmass(zR, aL ; n = 3)
    full(p::Float64) = Gmm(0.5,p) + aL*p*zratio*Gpm(0.5,p)
    pp = zeros(n)
    pp[1] = nsolve(full, 0 ; lower = 2., upper = 3.)
    if n > 1
        for i=2:n
            pp[i] = nsolve(full, 0 ; lower = pp[i-1]+2, upper = pp[i-1]+4)
        end
    end

    return pp./zR
end




# helper fuction for smooth
function flipdim(A::AbstractArray, d::Integer)
    nd = ndims(A)
    1 ≤ d ≤ nd || throw(ArgumentError("dimension $d is not 1 ≤ $d ≤ $nd"))
    if isempty(A)
        return copy(A)
    elseif nd == 1
        return reverse(A)
    end
    inds = indices(A)
    B = similar(A)
    nnd = 0
    for i = 1:nd
        nnd += Int(length(inds[i])==1 || i==d)
    end
    indsd = inds[d]
    sd = first(indsd)+last(indsd)
    if nnd==nd
        # flip along the only non-singleton dimension
        for i in indsd
            B[i] = A[sd-i]
        end
        return B
    end
    alli = [ indices(B,n) for n in 1:nd ]
    for i in indsd
        B[[ n==d ? sd-i : alli[n] for n in 1:nd ]...] = slicedim(A, d, i)
    end
    return B
end


function smooth(y ; win_len=11, win_method=2)
# This function requires the DSP package to be installed
# 1: flat
# 2: hanning
# 3: hamming ...
    if win_len%2==0
        win_len += 1 # only use odd numbers
    end
    if win_method == 1
        w = ones(win_len)
    elseif win_method==2
        w = DSP.hanning(win_len)
    elseif win_method==3
        w = DSP.hamming(win_len)
    end

    if win_len < 3
        return y
    elseif length(y) < win_len
        return y
    else
        y_new = [2*y[1].-flipdim(y[1:win_len],1); y[:]; 2*y[end].-flipdim(y[end-win_len:end],1)]
        y_smooth = conv(y_new, w/sum(w))
        ind = floor(Integer,1.5*win_len)
        return y_smooth[1+ind:end-ind-1]
    end

end;