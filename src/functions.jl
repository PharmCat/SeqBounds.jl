
#SeqBounds.jl
struct SeqBoundsResult
    sf::String
    v::Vector
    alpha::Float64
    zb::Vector
    asfv::Vector
    asfvdiff::Vector
    stdv::Vector
    sdproc::Vector
end

function obf(p, alpha)
    2(1 - cdf(Normal(), quantile(Normal(), 1-alpha/2)/sqrt(p)))
end

function bounds(v::Vector, alpha::Float64; h::Float64 = 0.05)
    if alpha ≥ 1 || alpha ≤ 0 error("alpha should be in range: 0 < alpha < 1!") end
    if v[end] > 1 error("Last value of v shoul be  ≤ 1!") end
    if length(v) > 1
        for i = 2:length(v)
            if !(v[i] > v[i-1])  error("v[i] shoul be > v[i-1]") end
        end
    end
    zninf       = -8.0
    side        = 1
    asfv        = obf.(v, alpha) # Alpha spending function values
    asfvdiff    = similar(asfv)
    vdiff       = similar(v)
    asfvdiff[1] = asfv[1]
    vdiff[1]    = v[1]
    for i = 2:length(v)
        vdiff[i]    = v[i] - v[i-1]
        asfvdiff[i] = asfv[i] - asfv[i-1]
    end
    stdv     = sqrt.(vdiff)
    sdproc   = sqrt.(v)
    za    = similar(asfv)
    zb    = similar(asfv)
    ya    = similar(asfv)
    yb    = similar(asfv)
    nints = zeros(Int, length(v))
    grids = Vector{StepRangeLen}(undef, length(v))
    #
    za[1] = zninf # For one side
    zb[1] = quantile(Normal(), 1-asfvdiff[1]/side)
    ya[1] = za[1]*stdv[1]
    yb[1] = zb[1]*stdv[1]

    nints[1] = ceil((yb[1] - ya[1])/(h*stdv[1]))
    grids[1] = range(ya[1], yb[1], length=nints[1] + 1)
    last     = pdf.(Normal(0, stdv[1]), grids[1])

    for i = 2:length(v)
        if i > 2
            # For next step
            nints[i-1] = ceil((yb[i-1]-ya[i-1])/(h*stdv[i-1]))
            grids[i-1] = range(ya[i-1], yb[i-1], length=nints[i-1]+1)
            last     = highordgrid(last, grids[i-1], grids[i-2], stdv[i-1])
        end

        za[i]    = zninf
        ya[i]    = zninf*sdproc[i]
        zerof    = x -> asfvdiff[i] - integrategrid(x, last, grids[i-1], stdv[i])

        yb[i]    = find_zero(zerof, (0.0, zb[i-1]))
        zb[i]    = yb[i]/sdproc[i]
    end
    SeqBoundsResult("O'Brien-Fleming",
    v,
    alpha,
    zb,
    asfv,
    asfvdiff,
    stdv,
    sdproc)
end

function integrategrid(x, last::Vector{T}, grid, stdv) where T
    hlast = step(grid)
    DIST  = Normal(x, stdv)
    fst   = cdf(DIST, grid[1]) * last[1]
    lst   = cdf(DIST, grid[end]) * last[end]
    sfn   = zero(T)
    for i = 2:length(grid)-1
        sfn += cdf(DIST, grid[i]) * last[i]
    end
    hlast*(sfn + fst/2 + lst/2)
end

function highordgrid(last::Vector{T}, x, grid, stdv) where T
    h    = step(grid)
    rvec = zeros(T, length(x))
    for n = 1:length(x)
        fst = last[1] * pdf(Normal(x[n], stdv), grid[1])
        lst = last[end] * pdf(Normal(x[n], stdv), grid[end])
        for m = 2:length(grid)-1
            rvec[n] += last[m] * pdf(Normal(x[n], stdv), grid[m])
        end
        rvec[n] += fst/2 + lst/2
    end
    rvec * h
end

function Base.show(io::IO, obj::SeqBoundsResult)
    println(io, """One-sided group sequential design
    Alpha spending function: $(obj.sf),  Alpha = $(obj.alpha)""")
    pretty_table(io, hcat(obj.v, obj.asfv, obj.asfvdiff, obj.zb, ccdf.(Normal(), obj.zb)); header = ["Portion", "Function value", "Spend", "Z", "Nominal p"])
end
################################################################################
#=
Dounds can be found by hcubature with HCubature.jl but it takes more time and memory
=#
#=
f1i(x) = quadgk(ab -> pdf(Normal(0, stdv[1]), ab), ya[1], yb[1], rtol=1e-8)[1]
f1i(yb[1])
f1(x) =  pdf(Normal(0, stdv[1]), x)

f21(x) = quadgk(u -> cdf(Normal(x, stdv[2]), u)*pdf(Normal(0, stdv[1]), u), ya[1], yb[1], rtol=1e-8)[1]
f21(yb[2])

# OR

f22(x) = quadgk(u -> cdf(Normal(x, stdv[2]), u)*f1(u), ya[1], yb[1], rtol=1e-8)[1]
zerof = x ->asfvdiff[2] - f23(x)
fz = find_zero(zerof, (0.0, zb[1]))
zbf = fz/sdproc[2]

f3(x) = hcubature(u -> cdf(Normal(x, stdv[i]), u[1]) * pdf(Normal(u[1], stdv[i-1]), u[2]) * pdf(Normal(0, stdv[1]), u[2]) , [ya[i-1], ya[i-2]], [yb[i-1],yb[i-2]], rtol=1e-8)[1]
f3(yb[3])
zerof = x ->asfvdiff[3] - f3(x)
fz = find_zero(zerof, (0.0, zb[2]))
zbf = fz/sdproc[3]
=#
