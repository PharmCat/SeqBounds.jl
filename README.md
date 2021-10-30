# SeqBounds

Group sequential design bounds.

## Install

```
import Pkg; Pkg.add("SeqBounds")
```

## Using

```
using SeqBounds

bounds([0.25, 0.5, 0.75, 1.0], 0.05; h = 0.05)
```

Result:

```
julia> bounds([0.25, 0.5, 0.75, 1.0], 0.05; h = 0.05)
One-sided group sequential design
Alpha spending function: O'Brien-Fleming,  Alpha = 0.05
┌─────────┬────────────────┬────────────┬─────────┬────────────┐
│ Portion │ Function value │      Spend │       Z │  Nominal p │
├─────────┼────────────────┼────────────┼─────────┼────────────┤
│    0.25 │     8.85754e-5 │ 8.85754e-5 │ 3.74955 │ 8.85754e-5 │
│     0.5 │      0.0055746 │ 0.00548602 │ 2.53993 │ 0.00554366 │
│    0.75 │      0.0236251 │  0.0180505 │ 2.01604 │  0.0218979 │
│     1.0 │           0.05 │  0.0263749 │ 1.72014 │  0.0427037 │
└─────────┴────────────────┴────────────┴─────────┴────────────┘
```

## API
```
  bounds(v::Vector, alpha::Float64; side = :one, asf = :obf, h::Float64 = 0.05)
```

Where:

* `v` - vector of information portion for each interim analysis;
* `alpha` - total alpha level;
* `h` - grid step multiplier, default 0.05, use 0.025 for better precision.

* `v` - vector of information portion for each interim analysis;
* `alpha` - total alpha level;
* `side` - one- or two- sided test (:one or :two);
* `asf` - alpha spending function (:obf, :pocock, :power, :ep);
* `h` - grid step multiplier, default 0.05, use 0.025 for better precision.

asf:

`:obf` - O'Brien-Fleming function

2(1 - cdf(Normal(), quantile(Normal(), 1 - α / side / 2)/sqrt(t)))

`:pocock` - Pocock function.

α / side * log(1 + (exp(1) - 1) * t)

`:power` - power function.

α / side * t^γ

`:ep` - equal parts.

## Reference

* Reboussin, D. M., DeMets, D. L., Kim, K., & Lan, K. K. G. (2000). Computations for Group Sequential Boundaries Using the Lan-DeMets Spending Function Method. Controlled Clinical Trials, 21(3), 190–207. doi:10.1016/s0197-2456(00)00057-x
* O’Brien, P. C., & Fleming, T. R. (1979). A Multiple Testing Procedure for Clinical Trials. Biometrics, 35(3), 549. doi:10.2307/2530245
* Lan KKG, DeMets DL. Discrete sequential boundaries for clinical trials. Biometrika. 1983; 70:659-63.
* DeMets DL, Lan G. “The alpha spending function approach to interim data analyses” in Recent Advances in Clinical Trial Design and Analysis, ed. PF Thall. Kluwer Academic Publishers, Boston; 1995.
* Armitage P, McPherson CK, Rowe BC. Repeated significance tests on accumulating data. Journal of the Royal Statistical Society. 1969; 132:235-44
* ldbounds: Lan-DeMets Method for Group Sequential Boundaries - https://CRAN.R-project.org/package=ldbounds (comparation)
