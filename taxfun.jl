function tau(y, state, params)
    eitcparams = params.eitcparams
    @unpack cutoff1, cutoff2, cutoff3, rate1, rate2 = eitcparams
    @unpack haskid = state
    if haskid == false
        return 0.
    end
    if y < 0
        return 0.
    elseif y >= 0 && y < cutoff1
        return rate1 * y
    elseif y >= cutoff1 && y < cutoff2
        return rate1 * cutoff1
    elseif y >= cutoff2 && y < cutoff3
        return rate1 * cutoff1 - rate2 * (y - cutoff2)
    else
        return 0.
    end
end
