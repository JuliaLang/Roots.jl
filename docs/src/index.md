# Roots.jl

Documentation for [Roots.jl](https://github.com/JuliaMath/Roots.jl)


## About

`Roots` is  a `Julia` package  for finding zeros of continuous
scalar functions of a single real variable. That  is solving $f(x)=0$ for $x$.
The `find_zero`function provides the
primary interface. It supports various algorithms through the
specification of a method. These include:

* Bisection-like algorithms. For functions where a bracketing interval
  is known (one where $f(a)$ and $f(b)$ have alternate signs), the
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
  
  
  
## Basic usage

Consider  the polynomial   function  $f(x) = x^5 - x + 1/2$. As a polynomial,  its roots, or  zeros, could  be identified with the  `roots` function of  the `Polynomials` package. However, even  that function uses a numeric method to identify   the values, as no  solution with radicals is available. That is, even for polynomials, non-linear root finders are needed to solve $f(x)=0$. 

The `Roots` package provides a variety algorithms for this  task. In this overview, only the  default ones  are illustrated.

For  the function $f(x) = x^5 - x + 1/2$ a simple plot will show a zero  somewhere between $-1.2$ and $-1.0$ and two zeros near $0.6$. 

For the zero between two values at which the function changes sign, a
bracketing method is useful, as bracketing methods are guaranteed to
converge for continuous functions by the intermediate value
theorem. A bracketing algorithm will be used when the initial data is
passed as a tuple:

```jldoctest find_zero
julia> using Roots

julia> f(x) =  x^5 - x + 1/2
f (generic function with 1 method)

julia> find_zero(f, (-1.2,  -1))
-1.0983313019186334
```

The default algorithm is guaranteed to have an  answer nearly as accurate as is  possible  given the limitations of floating point  computations. 

For the zeros "near" a point,  a non-bracketing method is often used, as generally  the algorithms are more efficient and can be  used in cases where a zero does  not. Passing just  the initial point will dispatch to  such a method:

```jldoctest find_zero
julia> find_zero(f,  0.6)
0.550606579334135
```


This finds  the answer  to the left of the starting point. To get the other nearby zero, a starting point closer to the answer cana be used.  However,  an initial graph might convince one  that any of the upto 5 reaal  roots  will   occur between `-5`  and `5`.  The `find_zeros` function uses  heuristics and a few of the  algorithms to  identify  all zeros between the specified range. Here  we see  there  are 3:

```jldoctest find_zero
julia> find_zeros(f, -5,  5)
3-element Vector{Float64}:
 -1.0983313019186334
  0.550606579334135
  0.7690997031778959
```
