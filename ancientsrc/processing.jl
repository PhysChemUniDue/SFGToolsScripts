function arr2mat(data::Array{SFSpectrum,1})
    hpix, vpix, nspec = size(data[1])
    n = length(data)
    m = zeros(n, hpix, vpix, nspec)
    faulty_spectra = Int64[]
    for i = 1:n, j=1:hpix, k=1:vpix, l=1:nspec
        # Check if that spectrum really was recorded
        if size(data[i], 3) < nspec
            push!(faulty_spectra, i)
            continue
        end
        m[i,j,k,l] = data[i][j,k,l]
    end
    unique!(faulty_spectra)
    for f in faulty_spectra
        @warn "faulty spectrum detected at index $f"
    end
    m
end
