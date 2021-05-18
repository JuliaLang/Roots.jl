# Attempt to find all zeros in an interval (a,b)

# Algorithm due to @djsegal in https://github.com/JuliaMath/Roots.jl/pull/113


# A naive approach to find zeros: split (a,b) by n points, look into each for a zero
# * k is oversampling rate for bisection. (It is relatively cheap to check for a bracket so we
#   oversample our intervals looking for brackets
# * assumes f(a) *not* a zero
function _fz(f, a, b, no_pts, k=4)
    zs = Real[]
    _fz!(zs, f, a, b, no_pts, k)
    zs
end

function _fz!(zs, f, a::T, b, no_pts, k=4) where {T}

    pts = range(a, stop=b, length=(no_pts-1)*k+1)
    n::Int = length(pts)

    fs = f.(pts)
    sfs = sign.(fs)

    u::T = first(pts)  # keep track of bigger interval
    found_bisection_zero = false

    for (i,x) in enumerate(pts[1:end])
        q,r = divrem(i-1, k)

        if i > 1 && iszero(r)
            v::T = x
            if !found_bisection_zero
                try
                    p1::T = identify_starting_point(u, v, sfs[(i-k):i])
                    rt::T = dfree(f, p1)
                    if !isnan(rt) && u < rt <= v
                        push!(zs, rt)
                    end
                catch err
                end
            end
            u = v
            found_bisection_zero = false
        end

        if i < n
            if iszero(fs[i+1])
                found_bisection_zero = true # kinda
                push!(zs, pts[i+1])
            elseif  sfs[i] * sfs[i+1] < 0
                found_bisection_zero = true
                rt = bisection(f, x, pts[i+1])
                push!(zs, rt)
            end
        end
    end

    sort!(zs)
end

# the algorithm first scans for zeros using the naive approach, then
# splits (a,b) by these zeros. This struct holds the subintervals
struct Interval{T}
a::T
b::T
depth::Int
end

Base.show(io::IO, alpha::Interval) = print(io, "($(alpha.a), $(alpha.b))")

# check if f(a) is non zero using tolerances max(atol, eps()), rtol
function _non_zero(fa, a::T, atol, rtol) where {T}
    abs(fa) >= max(atol, abs(a)*rtol*oneunit(fa)/oneunit(a), oneunit(fa)*eps(T))
end

# After splitting by zeros we have intervals (zm, zn) this is used to shrink
# to (zm+, zn-) where both are non-zeros, as defined above
function find_non_zero(f, a::T, barrier, xatol, xrtol, atol, rtol) where {T}
    nan = (0*a)/(0*a) # try to get typed NaN
    xtol = max(xatol, abs(a) * xrtol, oneunit(a) * eps(T))
    sgn = barrier > a ? 1 : -1
    ctr = 0
    x = a + 2^ctr*sgn*xtol
    while !_non_zero(f(x), x, atol, rtol)
        ctr += 1
        x += 2^ctr*sgn*xtol
        ((sgn > 0 && x > barrier) || (sgn < 0 && x < barrier)) && return nan
        ctr > 100 && return nan
    end

    x
end

# split a < z1 < z2 < ... < zn < b into intervals (a+,z1-), (z1+, z2-), ...
# where possible; push! onto ints
function make_intervals!(ints, f, a, b, zs, depth, xatol, xrtol, atol, rtol)
    pts = vcat(a, zs, b)

    for (u,v) in zip(pts[1:end-1], pts[2:end])
        ur = find_non_zero(f, u, v, xatol, xrtol, atol, rtol)
        isnan(ur) && continue

        vl = find_non_zero(f, v, u, xatol, xrtol, atol, rtol)
        isnan(vl) && continue

        push!(ints, Interval(ur, vl, depth))
    end
end


# adjust what we mean by x1 ~ x2 for purposes of adding a new zero
function approx_close(z1, z2, xatol, xrtol)
    tol = max(sqrt(xatol), max(abs(z1), abs(z2)) * sqrt(xrtol))
    abs(z1 - z2) < tol
end

# is proposed not near xs? (false and we add proposed)
function not_near(proposed, xs, xatol, xrtol)
    n = length(xs)
    n <= 1 && return true
    ind = n + 1
    for (i,rt) in enumerate(xs)
        if proposed < rt
            ind = i
            break
        end
    end
    if ind > 1 # xs[ind-1] <= propose < xs[ind]
        rt = xs[ind-1]
        approx_close(proposed, rt, xatol, xrtol) && return false
    end
    if ind <= n # value to right
        rt = xs[ind]
        approx_close(proposed, rt, xatol, xrtol) && return false
    end
    return true
end


"""
    find_zeros(f, a, b; [no_pts=12, k=8, naive=false, xatol, xrtol, atol, rtol])

Search for zeros of `f` in the interval `[a,b]`.


# Examples

```jldoctest find_zeros
julia> using Roots

julia> find_zeros(x -> exp(x) - x^4, -5, 20)        # a few well-spaced zeros
3-element Vector{Float64}:
 -0.8155534188089606
  1.4296118247255556
  8.613169456441398

julia> find_zeros(x -> sin(x^2) + cos(x)^2, 0, 2pi)  # many zeros
12-element Vector{Float64}:
 1.78518032659534
 2.391345462376604
 3.2852368649448853
 3.3625557095737544
 4.016412952618305
 4.325091924521049
 4.68952781386834
 5.00494459113514
 5.35145266881871
 5.552319796014526
 5.974560835055425
 6.039177477770888

julia> find_zeros(x -> cos(x) + cos(2x), 0, 4pi)    # mix of simple, non-simple zeros
6-element Vector{Float64}:
  1.0471975511965976
  3.141592653589943
  5.235987755982988
  7.330382858376184
  9.424777960769228
 11.519173063162574

julia> f(x) = (x-0.5) * (x-0.5001) * (x-1)          # nearby zeros
f (generic function with 1 method)

julia> find_zeros(f, 0, 2)
3-element Vector{Float64}:
 0.5
 0.5001
 1.0

julia> f(x) = (x-0.5) * (x-0.5001) * (x-4) * (x-4.001) * (x-4.2)
f (generic function with 1 method)

julia> find_zeros(f, 0, 10)
3-element Vector{Float64}:
 0.5
 0.5001
 4.2

julia> f(x) = (x-0.5)^2 * (x-0.5001)^3 * (x-4) * (x-4.001) * (x-4.2)^2  # hard to identify
f (generic function with 1 method)

julia> find_zeros(f, 0, 10, no_pts=21)                # too hard for default
5-element Vector{Float64}:
 0.49999999999999994
 0.5001
 4.0
 4.001
 4.200000000000001
```

!!! note
    There are two cases where the number of zeros may be underreported:
    * if the initial interval, [a,b], is too wide
    * if there are zeros  that are very nearby

----

The basic algorithm checks for zeros among the endpoints, and then
divides the interval `(a,b)` into `no_pts-1` subintervals and then
proceeds to look for zeros through bisection or a derivative-free
method.  As checking for a bracketing interval is relatively cheap and
bisection is guaranteed to converge, each interval has `k` pairs of
intermediate points checked for a bracket.

If any zeros are found, the algorithm uses these to partition `(a,b)`
into subintervals. Each subinterval is shrunk so that the endpoints
are not zeros and the process is repeated on the subinterval. If the
initial interval is too large, then the naive scanning for zeros may
be fruitless and no zeros will be reported. If there are nearby zeros,
the shrinking of the interval may jump over them, though as seen in
the examples, nearby roots can be identified correctly, though for
really nearby points, or very flat functions, it may help to increase
`no_pts`.


The tolerances are used to shrink the intervals, but not to find zeros
within a search. For searches, bisection is guaranteed to converge
with no specified tolerance. For the derivative free search, a
modification of the `Order0` method is used, which at worst case
compares `|f(x)| <= 8*eps(x)` to identify a zero. The algorithm might
identify more than one value for a zero, due to floating point
approximations. If a potential pair of zeros satisfy
`isapprox(a,b,atol=sqrt(xatol), rtol=sqrt(xrtol))` then they are
consolidated.

The algorithm can make many function calls. When zeros are found in an
interval, the naive search is carried out on each subinterval. To cut
down on function calls, though with some increased chance of missing
some zeros, the adaptive nature can be skipped with the argument
`naive=true` or the number of points stepped down.


The algorithm is derived from one in a
[PR](https://github.com/JuliaMath/Roots.jl/pull/113) by @djsegal.


!!! Note
The `IntervalRootFinding` package provides a rigorous alternative to this heuristic one. 
That package use interval arithmetic, so can compute bounds on the size of the image of
an interval under a function `f`. If this image includes `0`, then it can look for the zero. 
Bisection, on the other hand, only will look for a zero if the two endpoints have different signs,
a much more rigid condition for a potential zero.

For example, this function (due to `@truculentmath`) is particularly tricky, as it is positive at every floating point number, but has two zeros (the vertical asymptote at `15//11` is negative between adjacent floating point values):

```jldoctest find_zeros
julia> g(x) = x^2 + 1 +log(abs( 11*x-15 ))/99
g (generic function with 1 method)

julia> find_zeros(g, -3, 3)
Float64[]
```

The rigorous approach provided by `IntervalRootFinding` also has issues with this function, but it provides a more informed answer that the specified interval has an unknown number of zeros:


```
julia> using IntervalArithmetic, IntervalRootFinding

julia> g(x) = x^2 + 1 +log(abs( 11*x-15 ))/99
g (generic function with 1 method)

julia> IntervalRootFinding.roots(g, -3..3, IntervalRootFinding.Bisection)
1-element Vector{Root{Interval{Float64}}}:
 Root([1.36363, 1.36364], :unknown)
```

A less extreme usage might be the following, where `unique` indicates Bisection could be useful and indeed `find_zeros` will identify these values:

```
julia> g(x) = exp(x) - x^5
g (generic function with 1 method)

julia> rts = IntervalRootFinding.roots(g, -20..20)
2-element Vector{Root{Interval{Float64}}}:
 Root([12.7132, 12.7133], :unique)
 Root([1.29585, 1.29586], :unique)

julia> find_zeros(g, -20, 20)
2-element Vector{Float64}:
  1.2958555090953687
 12.713206788867632
```

"""
function find_zeros(f, a, b; no_pts = 12, k=8,
                    naive = false,
                    kwargs...
                    )

    a0, b0 = promote(float(a), float(b))
    a0 = isinf(a0) ? nextfloat(a0) : a0
    b0 = isinf(b0) ? prevfloat(b0) : b0

    # set tolerances if not specified
    fa0, fb0 = f(a0), f(b0)
    d = Dict(kwargs)
    T, S = eltype(a0), eltype(fa0)
    xatol::T = get(d, :xatol, eps(one(T))^(4//5) * oneunit(T))
    xrtol = get(d, :xrtol, eps(one(T)) * one(T))
    atol::S  = get(d, :atol,  eps(float(S)) * oneunit(S))
    rtol  = get(d, :rtol,  eps(float(S)) * one(S))

    zs = T[]  # collect zeros

    # check endpoints for exact zeros, then narrow
    abs(fa0) * oneunit(T) / oneunit(S) <= 8 * eps(a0) && push!(zs, a0)
    abs(fb0) * oneunit(T) / oneunit(S) <= 8 * eps(b0) && push!(zs, b0)
    a0 = find_non_zero(f, a0, b0, xatol, xrtol, atol, rtol)
    b0 = find_non_zero(f, b0, a0, xatol, xrtol, atol, rtol)


    _fz!(zs, f, a0, b0, no_pts,k)  # initial zeros

    ints = Interval{T}[] # collect subintervals
    !naive && !isempty(zs) &&  make_intervals!(ints, f, a0, b0, zs, 1, xatol, xrtol, atol, rtol)

    nzs = T[]
    cnt = 0

    while !naive && !isempty(ints)
        cnt += 1
        i = pop!(ints)
        # this is fussy. Ideally, we would want to explore intervals
        # with fewer points (each interval is already probed for k=4
        # bisections) but how many fewer?  We already had ~ (i.b-i.a) *
        # no_pts / (b-a) points. Would want an increased density so n
        # > (i.b - i.a) / (b - a) * no_pts but this doesn't perform as
        # well as we might expect

        # sub_no_pts = ceil(Int, (i.b - i.a) / (b-a) * no_pts * 2^(i.depth))
        #sub_no_pts <= 2 && continue  # stop on depth, always divide if roots
        #sub_no_pts = max(3, floor(Int, no_pts  / (2.0)^(i.depth)))
        sub_no_pts = floor(Int, no_pts  / (2.0)^(i.depth))

        empty!(nzs)
        if sub_no_pts >= 2
            _fz!(nzs, f, i.a, i.b, sub_no_pts, k)
        end

        if !isempty(nzs)
            azs = filter(rt->not_near(rt, zs, xatol, xrtol), nzs) # trim out nearby roots
            length(azs) == 0 && continue
            append!(zs, azs)
            sort!(zs)
            i.depth > 4 && continue
            make_intervals!(ints, f, i.a, i.b, azs, i.depth+1, xatol, xrtol, atol, rtol)
        end
    end

    length(zs) <= 1 && return zs
    sort!(zs)

    # may identify same zero with nearby values along the way
    # this trims out with a relaxed tolerance on how close
    # nearby roots can be. Default is epsilon^(2/5) ≈ 5e-7
    inds = Int[1]
    z1 = first(zs)
    for i in 2:length(zs)
        z2 = zs[i]
        if !approx_close(z1, z2, xatol, xrtol)
            push!(inds, i)
        end
        z1=z2
    end

    return zs[inds]

end
