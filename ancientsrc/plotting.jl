function heatmap(data::Array{SFSpectrum,1})
    m = arr2mat(data)
    ms = sum(m, dims=4)[:,:,:]
    a = []
    for i = 1:size(ms,3)
        f = figure()
        ax = subplot(111)
        title("Vertical Row $i")
        pcolormesh(ms[:,:,i])
        push!(a, ax)
    end
    a
end


function diffmap(data::Array{SFSpectrum,1}; refid=1, sigid=2)
    m = arr2mat(data)
    ms = sum(m, dims=4)[:,:,:]
    diffmat = ms[:,:,sigid] .- ms[:,:,refid]
    figure()
    ax = subplot(111)
    title("Vertical Row $sigid - $refid")
    pcolormesh(diffmat)
    ax
end
