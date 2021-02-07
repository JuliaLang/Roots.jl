[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliahub.com/docs/Roots/)
Linux: [![Build Status](https://travis-ci.org/JuliaMath/Roots.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Roots.jl)
Windows: [![Build status](https://ci.appveyor.com/api/projects/status/goteuptn5kypafyl?svg=true)](https://ci.appveyor.com/project/jverzani/roots-jl)

# Root finding functions for Julia


This package contains simple routines for finding roots of continuous
scalar functions of a single real variable. The `find_zero`function provides the
primary interface. It supports various algorithms through the
specification of a method. These include:

* Bisection-like algorithms. For functions where a bracketing interval
  is known (one where `f(a)` and `f(b)` have alternate signs), the
  `Bisection` method can be specified. For most floating point number
  types, bisection occurs in a manner exploiting floating point
  storage conventions. For others, an algorithm of Alefeld, Potra, and
  Shi is used. These methods are guaranteed to converge.


* Several derivative-free methods are implemented. These are specified
  through the methods `Order0`, `Order1` (the secant method), `Order2`
  (the Steffensen method), `Order5`, `Order8`, and `Order16`. The
  number indicates, roughly, the order of convergence. The `Order0`
  method is the default, and the most robust, but may take many more
  function calls to converge. The higher order methods promise higher
  order (faster) convergence, though don't always yield results with
  fewer function calls than `Order1` or `Order2`. The methods
  `Roots.Order1B` and `Roots.Order2B` are superlinear and quadratically converging
  methods independent of the multiplicity of the zero.


* There are historic methods that require a derivative or two:
  `Roots.Newton` and `Roots.Halley`.  `Roots.Schroder` provides a
  quadratic method, like Newton's method, which is independent of the
  multiplicity of the zero.

Each method's documentation has additional detail.

Some examples:


```julia
using Roots
f(x) = exp(x) - x^4

# a bisection method has the bracket specified with a tuple or vector
julia> find_zero(f, (8,9), Bisection())
8.613169456441398

julia> find_zero(f, (-10, 0))  # Bisection if x is a tuple and no method
-0.8155534188089606


julia> find_zero(f, (-10, 0), FalsePosition())  # just 11 function evaluations
-0.8155534188089607
```

For non-bracketing methods, the initial position is passed in as a
scalar:

```julia
## find_zero(f, x0::Number) will use Order0()
julia> find_zero(f, 3)         # default is Order0()
1.4296118247255556

julia> find_zero(f, 3, Order1()) # same answer, different method
1.4296118247255556

julia> find_zero(sin, BigFloat(3.0), Order16())
3.141592653589793238462643383279502884197169399375105820974944592307816406286198
```


The `find_zero` function can be used with callable objects:

```julia
using SymEngine
@vars x
find_zero(x^5 - x - 1, 1.0)  # 1.1673039782614185
```

Or,

```julia
using Polynomials
x = variable()
find_zero(x^5 - x - 1, 1.0)  # 1.1673039782614185
```

The function should respect the units of the `Unitful` package:

```julia
using Unitful
s = u"s"; m = u"m"
g = 9.8*m/s^2
v0 = 10m/s
y0 = 16m
y(t) = -g*t^2 + v0*t + y0
find_zero(y, 1s)      # 1.886053370668014 s
```

Newton's method can be used without taking derivatives, if the
`ForwardDiff` package is available:



```julia
using ForwardDiff
D(f) = x -> ForwardDiff.derivative(f,float(x))
```

Now we have:

```
f(x) = x^3 - 2x - 5
x0 = 2
find_zero((f,D(f)), x0, Roots.Newton())   # 2.0945514815423265
```

Automatic derivatives allow for easy solutions to finding critical
points of a function.

```julia
## mean
using Statistics
as = rand(5)

function M(x)
  sum([(x-a)^2 for a in as])
end

find_zero(D(M), .5) - mean(as)	  # 0.0

## median
function m(x)
  sum([abs(x-a) for a in as])
end

find_zero(D(m), (0, 1)) - median(as)	# 0.0
```

### Multiple zeros

The `find_zeros` function can be used to search for all zeros in a
specified interval. The basic algorithm essentially splits the interval into many
subintervals. For each, if there is a bracket, a bracketing algorithm
is used to identify a zero, otherwise a derivative free method is used
to search for zeros. This algorithm can miss zeros for various reasons, so the
results should be confirmed by other means.

```julia
f(x) = exp(x) - x^4
find_zeros(f, -10, 10)
```


### Convergence

For most algorithms, convergence is decided when

* The value |f(x_n)| < tol with `tol = max(atol, abs(x_n)*rtol)`, or

* the values x_n ≈ x_{n-1} with tolerances `xatol` and `xrtol` *and*
  f(x_n) ≈ 0 with a *relaxed* tolerance based on `atol` and `rtol`.

The algorithm stops if

* it encounters an `NaN` or an `Inf`, or

* the number of iterations exceed `maxevals`, or

* the number of function calls exceeds `maxfnevals`.

If the algorithm stops and the relaxed convergence criteria is met,
the suspected zero is returned. Otherwise an error is thrown
indicating no convergence. To adjust the tolerances, `find_zero`
accepts keyword arguments `atol`, `rtol`, `xatol`, and `xrtol`.


The `Bisection` and `Roots.A42` methods are guaranteed to converge
even if the tolerances are set to zero, so these are the
defaults. Non-zero values for `xatol` and `xrtol` can be specified to
reduce the number of function calls when lower precision is required.

Automatically `catch` convergence errors as follows:

```julia
f(x) = x ^ 2 + 2 # has no (real) root!

w = 0.0 # initialize root
try
    w = fzero(f, 0.001)
catch wuts_the_error
    if isa(wuts_the_error, Roots.ConvergenceFailed)
        w = NaN # assign root as NaN if convergence failed
    else
        error("something unrelated to convergence went wrong.")
    end
end
w # NaN
```


## An alternate interface

This functionality is provided by the `fzero` function, familiar to
MATLAB users. `Roots` also provides this alternative interface:


* `fzero(f, x0::Real; order=0)` calls a
  derivative-free method. with the order specifying one of
  `Order0`, `Order1`, etc.

* `fzero(f, a::Real, b::Real)` calls the `find_zero` algorithm with the
  `Bisection` method.

* `fzeros(f, a::Real, b::Real)` will call `find_zeros`.



## Usage examples

```julia
f(x) = exp(x) - x^4
## bracketing
fzero(f, 8, 9)		          # 8.613169456441398
fzero(f, -10, 0)		      # -0.8155534188089606
fzeros(f, -10, 10)            # -0.815553, 1.42961  and 8.61317

## use a derivative free method
fzero(f, 3)			          # 1.4296118247255558

## use a different order
fzero(sin, big(3), order=16)  # 3.141592653589793...
```

### Technical difference between find_zero and fzero

The `fzero` function is not identical to `find_zero`. When a function, `f`,
is passed to `find_zero` the code is specialized to the function `f`
which means the first use of `f` will be slower due to compilation,
but subsequent uses will be faster. For `fzero`, the code is not
specialized to the function `f`, so the story is reversed.
