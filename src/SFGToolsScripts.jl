module SFGToolsScripts

export load_spectra, pump_probe_process

include("userinterface.jl")

### userinterface

using TerminalMenus
using SFGTools
using Dates
using Statistics
using PyPlot

function get_spectra_choices(df, spectra_options; msg="Choose Spectra:")
    # here we use the default `pagesize` 10
    menu = MultiSelectMenu(spectra_options)

    # `request` returns a `Set` of selected indices
    # if the menu is canceled (ctrl-c or q), return an empty set
    spectra_choices = request(msg, menu) |> collect

    ids = df[:id][spectra_choices][1]
end

function get_spectra_options(date)
    df = list_spectra(date=yearmonthday(date), group=true)
    spectra_options = convert(Array{String,1}, df[:name])
    df, spectra_options
end

function get_date_choice(date_options)
    # `pagesize` is the number of items to be displayed at a time.
    #  The UI will scroll if the number of options is greater
    #   than the `pagesize`
    menu = RadioMenu(date_options)

    # `request` displays the menu and returns the index after the
    #   user has selected a choice
    choice = request("Choose Date:", menu)
end

function get_date_options()
    df = list_spectra()    
    dates = floor.(df[:date], Dates.Day)
    unique_dates = reverse(unique(dates), dims=1)
    date_options = string.(unique_dates)
    unique_dates, date_options
end

function SFGTools.load_spectra(;msg="Load Spectra:")
    unique_dates, date_options = get_date_options()
    date_choice = get_date_choice(date_options)
    date = unique_dates[date_choice]
    df, spectra_options = get_spectra_options(date)
    spectra_choices = get_spectra_choices(df, spectra_options; msg=msg) |> sort!
    data = load_spectra(spectra_choices)
end

function fieldcorrection!(spectrum; bias=[], dark=[], flat=[], darkflat=[])

    function rm_offset!(s, offset)
        offset_m = mean(offset, dims=3)
        for i = 1:size(s, 3)
            s[:,:,i] .-= offset_m[:,:,1]
        end
        s
    end

    function rm_offset(s, offset)
        offset_m = mean(offset, dims=3)
        n = deepcopy(s)
        for i = 1:size(s, 3)
            n[:,:,i] .-= offset_m[:,:,1]
        end
        n
    end

    function flatcorrection!(s, flat)
        flatmean = mean(flat, dims=3)
        flatmean ./= maximum(flatmean)
        for i = 1:size(s, 3)
            s[:,:,i] ./= flatmean[:,:,1]
        end
        s
    end

    # For readability make same length as all
    dafl = darkflat
    spec = spectrum

    !isempty(bias) &&            rm_offset!(spec, bias)
    !isempty(dark) && (dark_ub = rm_offset( dark, bias))
    !isempty(flat) && (flat_ub = rm_offset( flat, bias))
    !isempty(dafl) && (dafl_ub = rm_offset( dafl, bias))
    !isempty(dark) && rm_offset!(spec,    dark_ub)
    !isempty(dafl) && rm_offset!(flat_ub, dafl_ub)

    !isempty(flat) && flatcorrection!(spec, flat_ub)

    return spec
end

function fieldcorrectionreport(spectrum_raw, spectrum, bias, dark, flat, darkflat)
    pygui()
    λ = 1:size(spectrum[1], 1)#get_ir_wavelength(spectrum[1])
    figure()

    suptitle("Field Correction Report")

    subplot(3,2,1)
    pcolor(mean(spectrum_raw[end], dims=3)[:,:,1])
    colorbar()
    title("Last Spectrum")

    subplot(3,2,2)
    !isempty(bias) && pcolor(mean(bias, dims=3)[:,:,1])
    colorbar()
    title("Bias")

    subplot(3,2,3)
    !isempty(dark) && pcolor(mean(dark, dims=3)[:,:,1])
    colorbar()
    title("Dark")

    subplot(3,2,4)
    !isempty(flat) && pcolor(mean(flat, dims=3)[:,:,1])
    colorbar()
    title("Flat")

    subplot(3,2,5)
    !isempty(darkflat) && pcolor(mean(darkflat, dims=3)[:,:,1])
    colorbar()
    title("Dark Flat")

    subplot(3,2,6)
    plot(λ, mean(spectrum_raw[1], dims=2)[:,1,1], label="Raw")
    plot(λ, mean(spectrum[1],     dims=2)[:,1,1], label="Corrected")
    colorbar()
    title("Result")
    legend()

    plt[:tight_layout]()
    show()
    return nothing
end

function select_spectra()
    println("""Choose Field Correction Mode:
                  - none
               b - bias
               d - dark
               f - flat
               r - flat dark
               """)
    mode = readline()

    spectrum = load_spectra(msg="Load Spectra:")
    occursin("b", mode) ? (bias     = load_spectra(msg="Load Bias:")[1])      : (bias = [])
    occursin("d", mode) ? (dark     = load_spectra(msg="Load Dark:")[1])      : (dark = [])
    occursin("f", mode) ? (flat     = load_spectra(msg="Load Flat:")[1])      : (flat = [])
    occursin("r", mode) ? (darkflat = load_spectra(msg="Load Dark Flat:")[1]) : (darkflat = [])

    spectrum, bias, dark, flat, darkflat
end

function scanseries2mat(s::Array{SFSpectrum}, page=1)
    T = typeof(first(s).s)
    x = size(first(s), 1)
    y = size(first(s), 2)
    M = T(undef, length(s), x, y)
    for i = 1:length(s)
        M[i,:,:] = s[i][:,:,page]
    end
    return M
end

function plotheatmaps(m, dltime, sigcol, refcol)
    figure()

    x = 1:size(m,2)

    subplot(1,3,1)
    title("Signal")
    pcolor(x, dltime, m[:,:,sigcol])
    
    subplot(1,3,2)
    title("Reference")
    pcolor(x, dltime, m[:,:,refcol])

    subplot(1,3,3)
    title("Difference")
    pcolor(x, dltime, m[:,:,sigcol] .- m[:,:,refcol])

    return nothing
end

function plotintegrals(m, dltime, pixellocs, pixelwidth, sigcol, refcol)

    m_diff = m[:,:,sigcol] .- m[:,:,refcol]

    pixelstart = pixellocs .- pixelwidth ÷ 2
    pixelend   = pixellocs .+ pixelwidth ÷ 2
    
    figure()

    subplot(1,2,1)
    title("Spectrum")
    spec = mean(m[:,:,sigcol], dims=1)
    plot(1:size(spec, 2), spec[1,:])
    for i = 1:length(pixellocs)
        vlines([pixelstart[i], pixelend[i]], minimum(spec), maximum(spec), color="C$(i-1)")
    end


    subplot(1,2,2)
    title("Integrals")
    for i = 1:length(pixellocs)
        integ = mean(m_diff[:,pixelstart[i]:pixelend[i],1], dims=2)
        plot(dltime, integ[:,1,1])
    end
    legend()
    
    nothing
end

function pump_probe_process(sigcol=2, refcol=3)
    spectrum, bias, dark, flat, darkflat = select_spectra()

    # Remove Events
    rm_events!.(spectrum)
    !isempty(bias) && rm_events!(bias)
    !isempty(dark) && rm_events!(dark)
    !isempty(flat) && rm_events!(flat)
    !isempty(darkflat) && rm_events!(darkflat)

    # Field Correction
    spectrum_raw = deepcopy(spectrum)
    fieldcorrection!.(spectrum, bias=bias, dark=dark, flat=flat, darkflat=darkflat)
    fieldcorrectionreport(spectrum_raw, spectrum, bias, dark, flat, darkflat)

    # Average
    spectrum_avg = average.(spectrum)

    # Delay Time of Pump
    dltime = get_pump_delay.(spectrum)

    # Put spectrum in matrix
    m = scanseries2mat(spectrum)

    # Plot Heatmaps
    plotheatmaps(m, dltime, sigcol, refcol)

    # Plot Integrals
    println("Choose pixel locations over which to integrate.\nExample: [50, 80, 93]")
    pixellocs = readline() |> Meta.parse |> eval
    println("Choose width over which to integrate\nExample: 10")
    pixelwidth = readline() |> Meta.parse |> eval
    plotintegrals(m, dltime, pixellocs, pixelwidth, sigcol, refcol)
    
end

end # module
