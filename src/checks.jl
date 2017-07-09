""" Adds check for functional type """
macro check_functional(funcname, functype)
    msg = "Incorrect number of arguments: input is not an $functype functional"
    quote
        family($(esc(funcname))) ≠ Constants.$functype && throw(ArgumentError($msg))
    end
end

""" Adds argument check for energy derivative availability """
macro check_availability(funcname, functype)
    msg = "Functional does not implement energy."
    quote
        Constants.$functype ∉ flags($(esc(funcname))) && error($msg)
    end
end
""" Adds argument check for size compatibility """
macro check_size(funcname, rhoname, outname, factor)
    msg = "sizes of $rhoname and $outname are incompatible"
    quote
        if size($(esc(outname))) ≠ output_size($(esc(funcname)), $(esc(rhoname)), $factor)
            throw(ArgumentError($msg))
        end
    end
end

