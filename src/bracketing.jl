###


const bracketing_error = """The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
Consider a different bracket or try fzero(f, c) with an initial guess c.

"""

## Methods for root finding which use a bracket

## Bisection for FLoat64 values.
##
## From Jason Merrill https://gist.github.com/jwmerrill/9012954
## cf. http://squishythinking.com/2014/02/22/bisecting-floats/
# Alternative "mean" definition that operates on the binary representation
# of a float. Using this definition, bisection will never take more than
# 64 steps (over Float64)
const FloatNN = Union{Float64, Float32, Float16}



"""
    Roots.bisection64(f, a, b)

(unexported)

* `f`: a callable object, like a function

* `a`, `b`: Real values specifying a *bracketing* interval (one with
`f(a) * f(b) < 0`). These will be converted to `Float64` values.

Runs the bisection method using midpoints determined by a trick
leveraging 64-bit floating point numbers. After ensuring the
intermediate bracketing interval does not straddle 0, the "midpoint"
is half way between the two values onces converted to unsigned 64-bit
integers. This means no more than 64 steps will be taken, fewer if `a`
and `b` already share some bits.

The process is guaranteed to return a value `c` with `f(c)` one of
`0`, `Inf`, or `NaN`; *or* one of `f(prevfloat(c))*f(c) < 0` or
`f(c)*f(nextfloat(c)) > 0` holding. 

This function is a bit faster than the slightly more general 
`find_zero(f, [a,b], Bisection())` call.

Due to Jason Merrill.

"""
function bisection64(f, a::T, b::T) where {T <: Union{Float64, Float32, Float16}}
    u, v = promote(float(a), float(b))
    if v < u
        u,v = v,u
    end

    isinf(u) && (u = nextfloat(u))
    isinf(v) && (u = prevfloat(u))

    su, sv = sign(u), sign(v)
    
    if su * sv < 0
        # move to 0
        c = zero(u)
        sfu, sfc = sign(f(u)), sign(f(c))
        if sfu == sfc
            u =  c
        else
            v = c
        end
    end
    
    T == Float64 && return _bisection64(T, UInt64, f, u, v)
    T == Float32 && return _bisection64(T, UInt32, f, u, v)
    return _bisection64(T, UInt16, f, u, v)
end

## a,b same sign or zero, sfa * sfb < 0 is assumed
function _bisection64(T, S, f, a, b)
    nan = (0*a)/(0*a)
    negate = sign(a) < 0 ? -one(T) : one(T)

    ai, bi = reinterpret(S, abs(a)), reinterpret(S, abs(b))
    
    fa = f(a)
    iszero(fa) && return a
    sfa = sign(f(a))
    iszero(f(b)) && return b
    ai == bi && return nan
    
    mi = (ai + bi ) >> 1
    m = negate * reinterpret(T, mi)

    while a < m < b
        

        sfm = sign(f(m))
        iszero(sfm) && return m
        isnan(sfm) && return m

        if sfa * sfm < 0
            b, bi = m,  mi
        else
            a, sfa, ai = m, sfm, mi
        end

        mi = (ai + bi) >> 1
        m = negate * reinterpret(T, mi)
    end
    return m
end
        


####
## find_zero interface.
"""

    Bisection()

If possible, will use the bisection method over `Float64` values. The
bisection method starts with a bracketing interval `[a,b]` and splits
it into two intervals `[a,c]` and `[c,b]`, If `c` is not a zero, then
one of these two will be a bracketing interval and the process
continues. The computation of `c` is done by `_middle`, which
reinterprets floating point values as unsigned integers and splits
there. This method avoids floating point issues and when the
tolerances are set to zero (the default) guarantees a "best" solution
(one where a zero is found or the bracketing interval is of the type
`[a, nextfloat(a)]`).

When tolerances are given, this algorithm terminates when the midpoint
is approximately equal to an endpoint using absolute tolerance `xatol`
and relative tolerance `xrtol`. 

When a zero tolerance is given and the values are not `Float64`
values, this will call the `A42` method which has guaranteed convergence.
    
"""
struct Bisection <: AbstractBisection end  # either solvable or A42
struct BisectionExact <: AbstractBisection end

"""
    Roots.A42()

Bracketing method which finds the root of a continuous function within
a provided interval [a, b], without requiring derivatives. It is based
on algorithm 4.2 described in: 1. G. E. Alefeld, F. A. Potra, and
Y. Shi, "Algorithm 748: enclosing zeros of continuous functions," ACM
Trans. Math. Softw. 21, 327–344 (1995), DOI: 10.1145/210089.210111 .
[link](http://www.ams.org/journals/mcom/1993-61-204/S0025-5718-1993-1192965-2/S0025-5718-1993-1192965-2.pdf)


The default tolerances are: `xatol=zero(T)`, `xrtol=eps(T)`,
`maxevals=15`, and `maxfnevals=typemax(Int)`. There is no check made
using `atol`, `rtol` (for convergence in the `f(x)` direction.

"""
mutable struct A42 <: AbstractBisection end

## tracks for bisection, different, we show bracketing interval
function log_step(l::Tracks, M::AbstractBisection, state)
    push!(l.xs, state.xn0)
    push!(l.xs, state.xn1) # we store [ai,bi, ai+1, bi+1, ...]
end
function show_tracks(l::Tracks, M::AbstractBisection)
    xs = l.xs
    n = length(xs)
    for (i,j) in enumerate(1:2:(n-1))
        println(@sprintf("(%s, %s) = (% 18.16f, % 18.16f)", "a_$(i-1)", "b_$(i-1)", xs[j], xs[j+1]))
    end
    println("")
end
        
    

## helper function
function adjust_bracket(x0)
    u, v = float.(promote(x0...))
    if u > v
        u,v = v,u
    end
    isinf(u) && (u = nextfloat(u))
    isinf(v) && (v = prevfloat(v))
    u, v
end

function init_state(method::AbstractBisection, fs, x)
    length(x) > 1 || throw(ArgumentError(bracketing_error))

    x0, x1 = adjust_bracket(x)
    m = _middle(x0, x1)

    y0, y1, fm = sign.(promote(fs(x0), fs(x1), fs(m)))
    y0 * y1 > 0 && throw(ArgumentError("bracketing_error"))
    
    state = UnivariateZeroState(x1, x0, m,
                                y1, y0, fm,
                                0, 3,
                                false, false, false, false,
                                "")
    state

end

function init_state!(state::UnivariateZeroState{T,S}, ::AbstractBisection, fs, x::Union{Tuple, Vector}) where {T, S}
    x0, x1 = adjust_bracket(x)
    m::T = _middle(x0, x1)
    fx0::S, fx1::S, fm::S = sign(fs(x0)), sign(fs(x1)), sign(fs(m))
    init_state!(state, x1, x0, m, fx1, fx0, fm)
end

# for Bisection, the defaults are zero tolerances and strict=true
function init_options(::M,
                      state::UnivariateZeroState{T,S};
                      xatol=missing,
                      xrtol=missing,
                      atol=missing,
                      rtol=missing,
                      maxevals::Int=typemax(Int),
                      maxfnevals::Int=typemax(Int)) where {M <: Union{Bisection, BisectionExact,  A42}, T, S}

    ## Where we set defaults
    x1 = real(oneunit(state.xn1))
    fx1 = real(oneunit(float(state.fxn1)))
    strict = true

    # all are 0 by default
    options = UnivariateZeroOptions(ismissing(xatol) ? zero(x1) : xatol,       # unit of x
                                    ismissing(xrtol) ? zero(x1/oneunit(x1)) : xrtol,               # unitless
                                    ismissing(atol)  ? zero(fx1) : atol,  # units of f(x)
                                    ismissing(rtol)  ? zero(fx1/oneunit(fx1)) : rtol,            # unitless
                                    maxevals, maxfnevals, strict)

    options
end

function init_options!(options::UnivariateZeroOptions{Q,R,S,T}, ::Bisection) where {Q, R, S, T}
    options.xabstol = zero(Q)
    options.xreltol = zero(R)
    options.abstol = zero(S)
    options.reltol = zero(T)
    options.maxevals = typemax(Int)
    options.strict = true
end

## This uses _middle bisection Find zero using modified bisection
## method for FloatXX arguments.  This is guaranteed to take no more
## steps the bits of the type. The `a42` alternative usually has fewer
## iterations, but this seems to find the value with fewer function
## evaluations.
##
## This terminates when there is no more subdivision or function is zero

_middle(x::Float64, y::Float64) = _middle(Float64, UInt64, x, y)
_middle(x::Float32, y::Float32) = _middle(Float32, UInt32, x, y)
_middle(x::Float16, y::Float16) = _middle(Float16, UInt16, x, y)
_middle(x::Number, y::Number) = 0.5*x + 0.5 * y # fall back or non Floats

function _middle(T, S, x, y)
    # Use the usual float rules for combining non-finite numbers
    if !isfinite(x) || !isfinite(y)
        return x + y
    end
    # Always return 0.0 when inputs have opposite sign
    if sign(x) != sign(y) && !iszero(x) && !iszero(y)
        return zero(T)
    end
  
    negate = sign(x) < 0 || sign(y) < 0

    # do division over unsigned integers with bit shift
    xint = reinterpret(S, abs(x))
    yint = reinterpret(S, abs(y))
    mid = (xint + yint) >> 1

    # reinterpret in original floating point
    unsigned = reinterpret(T, mid)

    negate ? -unsigned : unsigned
end

function update_state(method::Union{Bisection,BisectionExact}, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T<:Number,S<:Number}


    y0 = o.fxn0
    m::T = o.m  
    ym::S = o.fm #sign(fs(m))
    incfn(o)

    if iszero(ym)
        o.message = "Exact zero found"
        o.xn1 = m
        o.fxn1= m
        o.x_converged = true
        return nothing
    end
    
    if y0 * ym < 0
        o.xn1, o.fxn1 = m, ym
    else
        o.xn0, o.fxn0 = m, ym
    end

    o.m = _middle(o.xn0, o.xn1)
    o.fm = sign(fs(o.m))
    return nothing

end

## convergence is much different here
## the method converges,
## as we bound between x0, nextfloat(x0) is not measured by eps(), but eps(x0)
function assess_convergence(method::Union{Bisection}, state::UnivariateZeroState{T,S}, options) where {T, S}

   
    state.x_converged && return true

    x0, x1, m::T = state.xn0, state.xn1, state.m

    if !(x0 < m < x1)
        state.x_converged = true
        return true
    end

    tol = max(options.xabstol, max(abs(x0), abs(x1)) * options.xreltol)
    if x1 - x0 > tol 
        return false
    end
    
    
    state.message = ""
    state.x_converged = true
    return true
end

# for exact convergence, we can skip some steps
function assess_convergence(method::BisectionExact, state::UnivariateZeroState{T,S}, options) where {T, S}

    state.x_converged && return true
    
    x0, m::T, x1 = state.xn0, state.m, state.xn1

    x0 < m < x1 && return false

    state.x_converged = true
    return true
end



## Bisection has special cases
## for FloatNN types, we have a slightly faster `bisection64` method
## for zero tolerance, we have either BisectionExact or A42 methods
## for non-zero tolerances, we have either a general Bisection or an A42
function find_zero(fs, x0, method::M;
                   tracks = NullTracks(),
                   verbose=false,
                   kwargs...) where {M <: Union{Bisection}}
    
    x = adjust_bracket(x0)
    T = eltype(x[1])
    F = callable_function(fs)
    state = init_state(method, F, x)
    options = init_options(method, state; kwargs...)
    tol = max(options.xabstol, maximum(abs.(x)) * options.xreltol)

    l = (verbose && isa(tracks, NullTracks)) ? Tracks(eltype(state.xn1)[], eltype(state.fxn1)[]) : tracks
    
    if iszero(tol)
        if T <: FloatNN
            !verbose && return bisection64(F, state.xn0, state.xn1) # speedier
            find_zero(BisectionExact(), F, options, state, l)
        else
            return find_zero(F, x, A42())
        end
    else
        find_zero(method, F, options, state, l)
    end

    if verbose
        show_trace(method, state, l)
    end
    
    state.xn1
    
end


###################################################
#
## A42
#
# Finds the root of a continuous function within a provided
# interval [a, b], without requiring derivatives. It is based on algorithm 4.2
# described in: 1. G. E. Alefeld, F. A. Potra, and Y. Shi, "Algorithm 748:
# enclosing zeros of continuous functions," ACM Trans. Math. Softw. 21,
# 327–344 (1995).
#
# Originally by John Travers

## put in utils?
@inline isbracket(fa,fb) = sign(fa) * sign(fb) < 0

# f[b,a]
@inline f_ab(a,b,fa,fb) = (fb - fa) / (b-a)

# f[a,b,d]
@inline function f_abd(a,b,d,fa,fb,fd)
    fab, fbd = f_ab(a,b,fa,fb), f_ab(b,d,fb,fd)
    (fbd - fab)/(d-a)
end

# a bit better than a - fa/f_ab
@inline secant_step(a, b, fa, fb) =  a - fa * (b - a) / (fb - fa)


# of (a,fa), (b,fb) choose pair where |f| is smallest
@inline choose_smallest(a, b, fa, fb) = abs(fa) < abs(fb) ? (a,fa) : (b,fb)

# assume fc != 0
## return a1,b1,d with a < a1 <  < b1 < b, d not there
## 
@inline function bracket!(state) #a,b,c, fa, fb, fc)
    if isbracket(state.fxn0, state.fm)
        # switch b, c
        state.xn1, state.m = state.m, state.xn1
        state.fxn1, state.fm = state.fm, state.fxn1
    else
        # switch a, c
        state.xn0,  state.m = state.m, state.xn0
        state.fxn0, state.fm = state.fm, state.fxn0
    end
end

@inline function bracket(a,b,c, fa, fb, fc)

    if isbracket(fa, fc)
        # switch b,c
        return (a,c,b, fa, fc, fb)
    else
        # switch a,c
        return (c,b,a, fc, fb, fa)
    end
end

# return c in (a+delta, b-delta)
# adds part of `bracket` from paper with `delta`
function newton_quadratic(a::T, b, d, fa, fb, fd, k::Int, delta=zero(T)) where {T}

    A = f_abd(a,b,d,fa,fb,fd)
    r = isbracket(A,fa) ? b : a
    
    # use quadratic step; if that fails, use secant step; if that fails, bisection
    if !(isnan(A) || isinf(A)) || !iszero(A)
        B = f_ab(a,b,fa,fb)

        dr = zero(r)
        for i in 1:k
            Pr = fa + B * (r-a) +  A * (r-a)*(r-b)
            Prp = (B + A*(2r - a - b))
            r -= Pr / Prp
        end
        if a+2delta < r < b - 2delta
            return r
        end
    end

    # try secant step
    r =  secant_step(a, b, fa, fb)

    if a + 2delta < r < b - 2delta
        return r 
    end
    
    return _middle(a, b) # is in paper r + sgn * 2 * delta
    
end

# state
function init_state(M::A42, f, xs) 
    u, v = promote(float(xs[1]), float(xs[2]))
    if u > v
        u, v = v, u
    end
    fu, fv = promote(f(u), f(v))

    state = UnivariateZeroState(v, u, v, ## x1, x0, m 
                                fv, fu, fv,
                                0, 2,
                                false, false, false, false,
                                "")
    
    init_state!(state, M, f, (u,v), false)
    state
end

# secant step, then bracket for initial setup
function init_state!(state::UnivariateZeroState{T,S}, ::A42, f, xs::Union{Tuple, Vector}, compute_fx=true) where {T, S}

    if !compute_fx
        u, v = state.xn0, state.xn1
        fu, fv = state.fxn0, state.fxn1
    else
        u, v = promote(float(xs[1]), float(xs[2]))
        if u > v
            u, v = v, u
        end
        fu, fv = f(u), f(v)
        state.fnevals = 2
        isbracket(fu, fv) || throw(ArgumentError(bracketing_error))
    end

    c::T = secant_step(u, v, fu, fv)
    fc::S = f(c)
    incfn(state)

    init_state!(state, v, u, c, fv, fu, fc)
    bracket!(state)

    return nothing
end

# for A42, the defaults are reltol=eps(), atol=0; 20 evals and strict=true
# this *basically* follows the tol in the paper (2|u|*rtol + atol)
function init_options(::A42,
                      state::UnivariateZeroState{T,S};
                      xatol=missing,
                      xrtol=missing,
                      atol=missing,
                      rtol=missing,
                      maxevals::Int=15,
                      maxfnevals::Int=typemax(Int)) where {T,S}

    strict=true
    options = UnivariateZeroOptions(ismissing(xatol) ? zero(T) : xatol,       # unit of x
                                    ismissing(xrtol) ? eps(one(T)) : xrtol,   # unitless
                                    ismissing(atol)  ? zero(S) : atol,  # units of f(x)
                                    ismissing(rtol)  ? zero(one(S)) : rtol,   # unitless
                                    maxevals, maxfnevals, strict)

    options
end

function init_options!(options::UnivariateZeroOptions{Q,R,S,T}, ::A42) where {Q, R, S, T}
    options.xabstol = zero(Q)
    options.xreltol = zero(one(R))
    options.abstol = zero(S)
    options.reltol = eps(one(T))
    options.maxevals = 15
    options.strict = true
end

function assess_convergence(method::A42, state::UnivariateZeroState{T,S}, options) where {T,S}

    (state.stopped || state.x_converged || state.f_converged) && return true
    if state.steps > options.maxevals
        state.stopped = true
        state.message *= "Too many steps taken. "
        return true
    end

    if state.fnevals > options.maxfnevals
        state.stopped=true
        state.message *= "Too many function evaluations taken. "
        return true
    end

    for (x,fx) in ((state.xn0, state.fxn0), (state.xn1, state.fxn1), (state.m, state.fm))
        if iszero(fx)
            state.f_converged = true
            state.xn1=x
            state.fxn1=fx
            state.message *= "Exact zero found. "
            return true
        end
    end

    a,b = state.xn0, state.xn1
    tol = max(options.xabstol, max(abs(a),abs(b)) * options.xreltol)

    if abs(b-a) <= 2tol
        # pick smallest of a,b,m
        u::T, fu::S = choose_smallest(a,b, state.fxn0, state.fxn1)
        x, fx = choose_smallest(u, state.m, fu, state.fm)
        state.xn1 = x
        state.fxn1 = fx
        state.x_converged = true
        return true
    end

    return false
end
    
function check_zero(::A42, state, c, fc)
    if isnan(c)
        state.stopped = true
        state.xn1 = c
        state.message *= "NaN encountered. "
        return true
    elseif isinf(c)
        state.stopped = true
        state.xn1 = c
        state.message *= "Inf encountered. "
        return true
    elseif iszero(fc)
        state.f_converged=true
        state.message *= "Exact zero found. "
        state.xn1 = c
        state.fxn1 = fc
        return true
    end
    return false
end

## 3, maybe 4, functions calls per step
function update_state(M::A42, f, state::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}
    
    a::T,b::T,d::T = state.xn0, state.xn1, state.m
    fa::S,fb::S,fd::S = state.fxn0, state.fxn1, state.fm

    mu = 0.5
    an, bn = a, b
    lambda = 0.7
    tole = max(options.xabstol, max(abs(a),abs(b)) * options.xreltol)
    delta = lambda * tole
    
    c::T = newton_quadratic(a, b, d, fa, fb, fd, 2, delta)
    fc::S = f(c)
    incfn(state)
    check_zero(M, state, c, fc) && return nothing

    a,b,d,fa,fb,fd = bracket(a,b,c,fa,fb,fc)
    
    c = newton_quadratic(a,b,d,fa,fb,fd, 3, delta)
    fc = f(c)
    incfn(state)    
    check_zero(M, state, c, fc) && return nothing
    
    a, b, d, fa, fb, fd = bracket(a, b, c, fa, fb,fc)
    
    u::T, fu::S = choose_smallest(a, b, fa, fb)
    c = u - 2 * fu * (b - a) / (fb - fa)
    if abs(c - u) > 0.5 * (b - a)
        c = _middle(a, b) 
    end
    fc = f(c)
    incfn(state)    
    check_zero(M, state, c, fc) && return nothing

    ahat::T, bhat::T, dhat::T, fahat::S, fbhat::S, fdhat::S = bracket(a, b, c, fa, fb, fc)
    if bhat - ahat < mu * (b - a) 
        #a, b, d, fa, fb, fd = ahat, b, dhat, fahat, fb, fdhat
        a, b, d, fa, fb, fd = ahat, bhat, dhat, fahat, fbhat, fdhat
    else
        m::T = _middle(ahat, bhat)
        fm::S = f(m)
        incfn(state)
        a, b, d, fa, fb, fd = bracket(ahat, bhat, m, fahat, fbhat, fm)
    end
    state.xn0, state.xn1, state.m = a, b, d
    state.fxn0, state.fxn1, state.fm = fa, fb, fd

    return nothing
end




# original used this, when possible. Possible speed up
# # approximate zero of f using inverse cubic interpolation
# # if the new guess is outside [a, b] we use a quadratic step instead
# # based on algorithm on page 333 of [1]
# function ipzero(f, a, fa, b, fb, c, fc, d, fd)

#     Q11 = (c - d)*fc/(fd - fc)
#     Q21 = (b - c)*fb/(fc - fb)
#     Q31 = (a - b)*fa/(fb - fa)
#     D21 = (b - c)*fc/(fc - fb)
#     D31 = (a - b)*fb/(fb - fa)
#     Q22 = (D21 - Q11)*fb/(fd - fb)
#     Q32 = (D31 - Q21)*fa/(fc - fa)
#     D32 = (D31 - Q21)*fc/(fc - fa)
#     Q33 = (D32 - Q22)*fa/(fd - fa)
#     c = a + (Q31 + Q32 + Q33)
#     if (c <= a) || (c >= b)
#         return newton_quadratic(f, a, fa, b, fb, d, fd, 3)
#     end
#     return c, f(c)
# end


## ----------------------------

"""

    FalsePosition()

Use the [false
position](https://en.wikipedia.org/wiki/False_position_method) method
to find a zero for the function `f` within the bracketing interval
`[a,b]`.

The false position method is a modified bisection method, where the
midpoint between `[a_k, b_k]` is chosen to be the intersection point
of the secant line with the x axis, and not the average between the
two values.

To speed up convergence for concave functions, this algorithm
implements the 12 reduction factors of Galdino (*A family of regula
falsi root-finding methods*). These are specified by number, as in
`FalsePosition(2)` or by one of three names `FalsePosition(:pegasus)`,
`FalsePosition(:illinois)`, or `FalsePosition(:anderson_bjork)` (the
default). The default choice has generally better performance than the
others, though there are exceptions.

For some problems, the number of function calls can be greater than
for the `bisection64` method, but generally this algorithm will make
fewer function calls.

Examples
```
find_zero(x -> x^5 - x - 1, [-2, 2], FalsePosition())
```
"""
struct FalsePosition{R} <: AbstractBisection end
FalsePosition(x=:anderson_bjork) = FalsePosition{x}()

function update_state(method::FalsePosition, fs, o::UnivariateZeroState{T,S}, options::UnivariateZeroOptions) where {T,S}

    a::T, b::T =  o.xn0, o.xn1

    fa::S, fb::S = o.fxn0, o.fxn1

    lambda = fb / (fb - fa)
    tau = 1e-10                   # some engineering to avoid short moves
    if !(tau < abs(lambda) < 1-tau)
        lambda = 1/2
    end
    x::T = b - lambda * (b-a)        
    fx::S = fs(x)
    incfn(o)

    if iszero(fx)
        o.xn1 = x
        o.fxn1 = fx
        return
    end

    if sign(fx)*sign(fb) < 0
        a, fa = b, fb
    else
        fa = galdino_reduction(method, fa, fb, fx) #galdino[method.reduction_factor](fa, fb, fx)
    end
    b, fb = x, fx

    o.xn0, o.xn1 = a, b 
    o.fxn0, o.fxn1 = fa, fb
    
    nothing
end

# the 12 reduction factors offered by Galadino
galdino = Dict{Union{Int,Symbol},Function}(:1 => (fa, fb, fx) -> fa*fb/(fb+fx),
                                           :2 => (fa, fb, fx) -> (fa - fb)/2,
                                           :3 => (fa, fb, fx) -> (fa - fx)/(2 + fx/fb),
                                           :4 => (fa, fb, fx) -> (fa - fx)/(1 + fx/fb)^2,
                                           :5 => (fa, fb, fx) -> (fa -fx)/(1.5 + fx/fb)^2,
                                           :6 => (fa, fb, fx) -> (fa - fx)/(2 + fx/fb)^2,
                                           :7 => (fa, fb, fx) -> (fa + fx)/(2 + fx/fb)^2,
                                           :8 => (fa, fb, fx) -> fa/2,
                                           :9 => (fa, fb, fx) -> fa/(1 + fx/fb)^2,
                                           :10 => (fa, fb, fx) -> (fa-fx)/4,
                                           :11 => (fa, fb, fx) -> fx*fa/(fb+fx),
                                           :12 => (fa, fb, fx) -> (fa * (1-fx/fb > 0 ? 1-fx/fb : 1/2))  
)


# give common names
for (nm, i) in [(:pegasus, 1), (:illinois, 8), (:anderson_bjork, 12)]
    galdino[nm] = galdino[i]
end

# from Chris Elrod; https://raw.githubusercontent.com/chriselrod/AsymptoticPosteriors.jl/master/src/false_position.jl
@generated function galdino_reduction(methods::FalsePosition{R}, fa, fb, fx) where {R}
    f = galdino[R]
    quote
        $Expr(:meta, :inline)
        $f(fa, fb, fx)
    end
end
