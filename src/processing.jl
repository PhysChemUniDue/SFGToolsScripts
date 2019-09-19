using LsqFit, Statistics
using Dierckx

"""
`baseline_correction(spectrum; fitrange=(50, 50), poly=2)`
"""
function baseline_correction(spectrum; fitrange=(50, 50), poly=2)
    s = deepcopy(spectrum)

    function model(x,p)
        y = zeros(length(x))
        for n = 1:length(p)
            @. y += p[n] * x^(n-1)
        end
        y
    end

    n_rows = size(s,2)
    n_frames = size(s,3)
    n_pixels = size(s,1)
    x_data = [1:fitrange[1]..., n_pixels-fitrange[2]:n_pixels...]
    # Fit the baseline for each row and each frame individually
    baseline = zeros(n_pixels)
    for r in 1:n_rows, f in 1:n_frames
        y_data = s[x_data,r,f]
        p0 = zeros(poly)
        fit = curve_fit(model, x_data, y_data, p0)
        baseline .= model(1:n_pixels, fit.param)
        s[:,r,f] .-= baseline
    end
    s
end


"""
`get_reffactors(r; refrange=1:10, verbose=false, threshold=threshold=(0.2 * maximum(r)))`
"""
function get_reffactors(r; refrange=1:10, verbose=false, threshold=threshold=(0.2 * maximum(r)))
    b = r .> threshold
    m = fill(NaN, size(r)) # contains values above threshold
    m[b] .= r[b]

    mref = mean(r[refrange,:], dims=1)[:]  # reference spectrum
    mfac = similar(r) # factors how much the data point deviates from mref
    mfacfit = similar(mfac)

    for i in 1:size(r,1)
        mfac[i,:] .= m[i,:] ./ mref
    end

    x = 1:size(r,2)
    for i = 1:size(r,1)
        yf = mfac[i, .!isnan.(mfac[i,:])]
        xf = x[.!isnan.(mfac[i,:])]
        weights = m[i,:][.!isnan.(m[i,:])].^2
        weights ./= maximum(weights)
        spl = Spline1D(xf, yf, w=weights, s=5e-4, k=1, bc="nearest")
        mfacfit[i,:] .= evaluate(spl, x)
    end

    mfacfit
end

"""
normalize(r; navg=size(r,1))
"""
function normalize(r; navg=size(r,1))
    rs = sort(r[:])
    rmax_avg = mean(rs[end-navg+1:end])
    r ./ rmax_avg
end

"""
function subtract_diff(ref, sig; refrange=1:10)
"""
function subtract_diff(ref, sig; refrange=1:10)
    dif = sig .- ref
    mdiff = mean(dif[refrange,:], dims=1)[:]
    y = similar(sig)
    for i = 1:size(sig,1)
        y[i,:] .= sig[i,:] .- mdiff
    end
    y
end
