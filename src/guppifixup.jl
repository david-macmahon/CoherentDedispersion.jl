"""
    fixup!(grh::GuppiRaw.Header)

This function will fix broken/missing fields in `grh`, making some possibly
incorrect assumptions.  Header field `:name` may be excluded from being fixed up
by passing `fix_name=false` as a keyword argument.  Each field fixed by this
function will elicit a warning message.  The following fields have fixes defined
for them:

- `:obsbw` Set to -187.5 if missing
- `:obsbw` Forced to be negative if positive
- `:chan_bw` Set to OBSBW/(OBSNCHAN/NANTS) if missing
- `:chan_bw` Forced to be negative if positive
- `:ra` Converted from `"HH:MM:SS.s"` to numerical degrees if it is a `String`
- `:ra` Converted from `:ra_str` to numerical degrees if it is not present
- `:dec` Converted from `"DD:MM:SS.s"` to numerical degrees if it is a `String`
- `:dec` Converted from `:dec_str` to numerical degrees if it is not present
"""
function fixup!(grh::GuppiRaw.Header; kwargs...)
    if get(kwargs, :fix_obsbw, true)
        if !haskey(grh, :obsbw)
            @warn "setting obsbw to -187.5"
            grh[:obsbw] = -187.5
        elseif grh[:obsbw] > 0
            @warn "forcing obsbw to be negative"
            grh[:obsbw] = -grh[:obsbw]
        end
    end

    if get(kwargs, :fix_chan_bw, true)
        if !haskey(grh, :chan_bw)
            grh[:chan_bw] = grh[:obsbw] / (grh[:obsnchan]/get(grh,:nants,1))
        elseif grh[:chan_bw] > 0
            grh[:chan_bw] = -grh[:chan_bw]
        end
    end

    if get(kwargs, :fix_ra, true)
        if !haskey(grh, :ra)
            if !haskey(grh, :ra_str)
                @warn "setting ra to 0.0"
                grh[:ra] = 0.0
            else
                ra = hms2deg(grh[:ra_str])
                @warn "setting ra to $ra"
                grh[:ra] = ra
            end
        elseif grh[:ra] isa AbstractString
            ra = hms2deg(grh[:ra])
            @warn "setting ra to $ra"
            grh[:ra] = ra
        end
    end

    if get(kwargs, :fix_dec, true)
        if !haskey(grh, :dec)
            if !haskey(grh, :dec_str)
                @warn "setting dec to 0.0"
                grh[:dec] = 0.0
            else
                dec = dms2deg(grh[:dec_str])
                @warn "setting dec to $dec"
                grh[:dec] = dec
            end
        elseif grh[:dec] isa AbstractString
            dec = dms2deg(grh[:dec])
            @warn "setting dec to $dec"
            grh[:dec] = dec
        end
    end

    grh
end
