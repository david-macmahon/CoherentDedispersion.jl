# TODO should nintegrate come from consumed item?
function _outputtask(pqin; outbuf, nintegrate=4)
    outputfiles = Dict{String, IO}()

    # intvw is view into outbuf with extra dim for integration
    intvw = reshape(outbuf, 4, :, size(outbuf, 2))

    # intbuf is similar to outvw, but has singleton summation dimension
    intbuf = similar(outbuf, 1, size(outbuf,1)÷nintegrate, size(outbuf,2))

    # While non-empty items keep coming in
    while true
        consume!(pqin) do item::NamedTuple
            # If item is empty, return nothing
            isempty(item) && return nothing

            # Get dm, header, and block from item
            @show keys(item)
            dm = item.dm
            overlap = item.overlap
            rawname = item.rawname
            header = item.header
            data_blocks = item.data
            @info "outputtask got overlap $overlap pktidx $(
                header[:pktidx]) data $(
                size(data_blocks[1])) $(data_blocks[1][1])"

            # Compute Stokes I (for now)
            outbuf .= abs2.(data_blocks[1])
            outbuf .+= abs2.(data_blocks[2])

            # Integrate
            sum!(intbuf, intvw)
            nout = (size(outbuf, 1) - overlap) ÷ nintegrate

            # Construct outnane
            outname = replace(rawname,
                r"(\.\d\d\d\d)?\.raw$"=>".coddspec.$(lpad(nintegrate, 4, '0')).fil")

            # Get IO for outname
            io = get!(outputfiles, outname) do
                @info "opening $outname"
                newio = open(outname, "w")
                fbh = Filterbank.Header(header;
                    tsamp=nintegrate*header[:tbin], refdm=dm
                )
                write(newio, fbh)
                newio
            end

            # Write data
            write(io, @view intbuf[1, 1:nout, :])

            return item.data
        end === nothing && break
    end

    for (fn, io) in outputfiles
        @info "closing $fn"
        close(io)
    end

    @info "outputtask done"
    return nothing
end
