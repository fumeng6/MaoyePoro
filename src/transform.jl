# 源历史指流量、荷载等量随时间的变化，按真实时刻及 SI 单位解释；重复节点保留阶跃，窗外节点保留斜坡。
# 主求解使用 weighted_time_response，按每项加载开始后的延迟计算阶跃/斜坡响应。
# source_increments/time_response 是将线性变化离散为阶梯加载的辅助计算，采用有限时间步近似。
function source_increments(input::PoroelasticInput, dtg)
    dt = input.twindow/(input.nt-1)
    times, increments = Float64[], Float64[]
    if input.sinj0[1] != 0
        push!(times, 0.0)
        push!(increments, input.sinj0[1])
    end
    for i in 2:length(input.ts0)
        t0, t1 = input.ts0[i-1], input.ts0[i]
        dv = input.sinj0[i]-input.sinj0[i-1]
        dv == 0 && continue
        if t1 == t0
            push!(times, t1)
            push!(increments, dv)
            times[end] >= input.twindow && return times,increments
        else
            n = 1 + trunc(Int, (t1-t0)/min(dt, dtg))
            for j in 1:n
                push!(times, t0+j*(t1-t0)/n)
                push!(increments, dv/n)
                times[end] >= input.twindow && return times,increments
            end
        end
    end
    return times, increments
end

"""Right-continuous value of the supplied piecewise-linear source history."""
function source_value(input::PoroelasticInput, t)
    t < first(input.ts0) && return 0.0
    i = searchsortedlast(input.ts0,t)
    i == length(input.ts0) && return input.sinj0[end]
    a = (t-input.ts0[i])/(input.ts0[i+1]-input.ts0[i])
    return (1-a)*input.sinj0[i]+a*input.sinj0[i+1]
end

"""Invert a Laplace response and convolve its step response with a staircase history.

`direct` is the real high-frequency limit. The transform convention is
f(t)=∫F(f)exp(+2πift)df. `intdif=1` applies a node forward difference for
volume-change sources; `-1` integrates the response.
"""
function time_response(spectrum, dfg, alpha, input::PoroelasticInput, direct;
                       intdif=(input.istype == 1 ? 1 : 0), increments=nothing,
                       exact_direct=true)
    nfg = length(spectrum)
    ntg = 2nfg
    dtg = 1/(ntg*dfg)
    dt = input.twindow/(input.nt-1)
    c = real(direct)
    f = zeros(ComplexF64, ntg)
    # Dividing the residual transfer by s gives its step response in the
    # Laplace domain.
    for lf in 1:nfg
        s = complex(alpha,2pi*(lf-1)*dfg)
        f[lf] = (spectrum[lf]-c)/s
    end
    f[nfg+2:end] .= conj.(reverse(f[2:nfg]))
    step = real.(bfft(f)) .* dfg .* exp.(alpha .* ((0:ntg-1) .* dtg))
    # For this discrete Bromwich sum, the time integral from zero is
    # sum(R(s_l)/s_l * (exp(s_l*t)-1))*dfg. Subtracting the zero-time value
    # enforces that initial condition; frequency truncation and aliasing remain.
    step .-= step[1]
    ts, ds = isnothing(increments) ? source_increments(input, dtg) : increments
    out = zeros(Float64, input.nt)
    for it in eachindex(out)
        t = (it-1)*dt
        # c*δ(t) convolved with a linear history is exactly c*history(t).
        exact_direct && (out[it] = c*source_value(input,t))
        for j in eachindex(ts)
            ts[j] > t && continue
            x = (t-ts[j])/dtg
            i = min(1+floor(Int, x), ntg)
            b = x-floor(x)
            out[it] += ds[j]*((exact_direct ? 0.0 : c)+(1-b)*step[i]+b*step[min(i+1, ntg)])
        end
    end
    if intdif == 1
        if input.nt == 2
            out .= (out[2]-out[1])/dt
        elseif input.nt == 3
            d1, d2 = (out[2]-out[1])/dt, (out[3]-out[2])/dt
            out .= (d1,d2,d2)
        else
            for it in 2:input.nt-1
                out[it] = (out[it+1]-out[it])/dt
            end
            out[1] = 2out[2]-out[3]
            out[end] = out[end-1]
        end
    elseif intdif == -1
        out[1] = 0
        for it in 2:input.nt
            out[it] = out[it-1]+out[it]*dt
        end
    end
    return out
end

"""Exact Laplace transform of the supplied piecewise-linear history.

Repeated times retain their net jumps; after the final knot the value is held.
The transform requires nonzero s (the Bromwich line has Re(s)>0).
"""
function source_laplace(input::PoroelasticInput, s)
    iszero(s) && throw(ArgumentError("source_laplace requires nonzero s"))
    # Integrate each actual linear segment. Summing derivative jumps first
    # would cancel to O(s) for a finite pulse and lose precision near s=0.
    value = zero(s)
    for i in 2:length(input.ts0)
        t0,t1 = input.ts0[i-1],input.ts0[i]
        dt = t1-t0
        iszero(dt) && continue
        x = s*dt
        phi1 = iszero(x) ? one(x) : -expm1(-x)/x
        if abs(x) < 0.1
            term = one(x)/2
            phi2 = term
            for n in 1:8
                term *= -x*(n+1)/(n*(n+2))
                phi2 += term
            end
        else
            phi2 = (-expm1(-x)-x*exp(-x))/(x*x)
        end
        value += exp(-s*t0)*dt*(input.sinj0[i-1]*phi1+
                    (input.sinj0[i]-input.sinj0[i-1])*phi2)
    end
    value+input.sinj0[end]*exp(-s*input.ts0[end])/s
end

"""Net jumps of h, including the possible initial jump from zero.

If a transfer contains a linear*s term, these amounts times that coefficient
are Dirac impulse weights in the response, not finite node rates.
"""
function source_jump_amounts(input::PoroelasticInput)
    times,amounts = Float64[],Float64[]
    if input.sinj0[1] != 0
        push!(times,0.0); push!(amounts,input.sinj0[1])
    end
    for i in 2:length(input.ts0)
        input.ts0[i] == input.ts0[i-1] || continue
        delta = input.sinj0[i]-input.sinj0[i-1]
        iszero(delta) && continue
        if !isempty(times) && times[end] == input.ts0[i]
            amounts[end] += delta
        else
            push!(times,input.ts0[i]); push!(amounts,delta)
        end
    end
    [(time=t,amount=a) for (t,a) in zip(times,amounts) if a != 0]
end

"""Regular derivative of h: right limits internally, left limit at window end.

The Dirac weights of instantaneous jumps are returned by source_jump_amounts.
At a knot coinciding with the output-window end, use the slope of the segment
ending there; at interior knots, use the slope of the following segment.
"""
function source_slope(input::PoroelasticInput,t)
    t < first(input.ts0) && return 0.0
    i = searchsortedlast(input.ts0,t)
    if t == input.twindow && i > 0 && input.ts0[i] == t
        i = searchsortedfirst(input.ts0,t)-1
    end
    (i < 1 || i >= length(input.ts0)) && return 0.0
    (input.sinj0[i+1]-input.sinj0[i])/(input.ts0[i+1]-input.ts0[i])
end

"""Causal half derivative D^(1/2)h, including physical jump singularities."""
function source_half_derivative(input::PoroelasticInput,t)
    t < 0 && return 0.0
    current = source_value(input,t)
    first_at_t = searchsortedfirst(input.ts0,t)
    if first_at_t <= length(input.ts0) && input.ts0[first_at_t] == t
        left = first_at_t == 1 ? 0.0 : input.sinj0[first_at_t]
        jump = current-left
        iszero(jump) || return copysign(Inf,jump)
    end
    iszero(t) && return 0.0
    # Integrate the supplied history relative to its current value:
    # D^(1/2)h = [h(t)/sqrt(t) + 1/2∫(h(t)-h(u))/(t-u)^(3/2)du]/sqrt(pi).
    # This avoids cancellation of large ramp responses long after a short pulse.
    value = current/sqrt(t)
    for i in 2:length(input.ts0)
        t0,t1 = input.ts0[i-1],input.ts0[i]
        (t1 == t0 || t0 >= t) && continue
        finish = min(t1,t)
        a,b = sqrt(t-t0),sqrt(t-finish)
        factor = (finish-t0)/(a+b)^2
        value += factor*(current-input.sinj0[i-1])/a
        # At a continuous upper endpoint h(t)-h(t)=0; do not form 0/0.
        iszero(b) || (value += factor*(current-input.sinj0[i])/b)
    end
    value/sqrt(pi)
end

"""Causal half integral I^(1/2)h, evaluated directly on linear history pieces."""
function source_half_integral(input::PoroelasticInput,t)
    t < 0 && return 0.0
    value = 0.0
    for i in 2:length(input.ts0)
        t0,t1 = input.ts0[i-1],input.ts0[i]
        (t1 == t0 || t0 >= t) && continue
        finish = min(t1,t)
        a,b = sqrt(t-t0),sqrt(t-finish)
        h0 = input.sinj0[i-1]
        h1 = finish == t1 ? input.sinj0[i] :
            h0+(input.sinj0[i]-h0)*(finish-t0)/(t1-t0)
        # Exact endpoint weights for ∫h(u)/sqrt(t-u)du, expressed without
        # subtracting nearly equal roots or powers of the elapsed time.
        value += 2(finish-t0)/(3(a+b)^2)*((a+2b)*h0+(2a+b)*h1)
    end
    t > input.ts0[end] && (value += 2input.sinj0[end]*sqrt(t-input.ts0[end]))
    value/sqrt(pi)
end

"""Bound the causal history weight of one retained Fourier--Laplace sample.

The finite residual response uses step and ramp primitives starting at each
actual loading time. This triangle bound also covers their grid interpolation;
`dtg` is the largest interpolation step used by the refinement sequence.
Future-only history pieces have exactly zero weight.
"""
function causal_history_bound(input::PoroelasticInput,s,dtg)
    T = input.twindow
    a,z = real(s),abs(s)
    step(lag) = lag <= 0 ? 0.0 : (exp(a*(lag+dtg))+1)/z
    ramp(lag) = lag <= 0 ? 0.0 : abs(-expm1(-s*T))*(exp(a*(lag+dtg))+1)/z^2
    bound = abs(input.sinj0[1])*step(T)
    for j in 2:length(input.ts0)
        t0,t1 = input.ts0[j-1],input.ts0[j]
        dh = input.sinj0[j]-input.sinj0[j-1]
        iszero(dh) && continue
        bound += t1 == t0 ? abs(dh)*step(T-t1) :
            abs(dh/(t1-t0))*(ramp(T-t0)+ramp(T-t1))
    end
    bound
end

"""Invert a transfer and apply the history at its actual causal loading times.

G=linear*s+half*sqrt(s)+direct+inverse_half/sqrt(s)+R. The four elementary
terms act analytically on the exact history. The residual is integrated into
step and ramp primitives, each explicitly zero at nonpositive lag. Thus a
future load cannot enter through the period images of a global H(s) inversion.
The numerical ramp is held constant after lag T=input.twindow. Every requested
lag is at most T, so this is the same physical ramp on the needed interval;
the bounded continuation reduces positive-lag period-image contamination.
The residual step and ramp responses start at zero.

Returned samples contain the regular rate for the linear*s term. Its impulse
weights are linear*source_jump_amounts(input). A sqrt(s) term at a true history
jump returns the signed infinity of the causal half derivative.
"""
function weighted_time_response(Gs,dfg,alpha,input::PoroelasticInput;
        coeff=(linear=0.0,half=0.0,direct=0.0,inverse_half=0.0),zero_initial=true,
        fft_plan=nothing,fft_workspace=nothing)
    nfg = length(Gs)
    ntg = 2nfg
    dtg = 1/(ntg*dfg)
    linear = Float64(coeff.linear)
    half = Float64(coeff.half)
    direct = Float64(coeff.direct)
    inverse_half = hasproperty(coeff,:inverse_half) ? Float64(coeff.inverse_half) : 0.0
    needs_step = input.sinj0[1] != 0 || any(j ->
        input.ts0[j] == input.ts0[j-1] && input.ts0[j] < input.twindow &&
        input.sinj0[j] != input.sinj0[j-1], 2:length(input.ts0))
    needs_ramp = any(j -> input.ts0[j] > input.ts0[j-1] &&
        input.ts0[j-1] < input.twindow && input.sinj0[j] != input.sinj0[j-1],
        2:length(input.ts0))
    fs = needs_step ? (fft_workspace === nothing ? zeros(ComplexF64,ntg) :
        fill!(fft_workspace.fs,0)) : ComplexF64[]
    fa = needs_ramp ? (fft_workspace === nothing ? zeros(ComplexF64,ntg) :
        fill!(fft_workspace.fa,0)) : ComplexF64[]
    for lf in 1:nfg
        s = complex(alpha,2pi*(lf-1)*dfg)
        residual = Gs[lf]-linear*s-half*sqrt(s)-direct-inverse_half/sqrt(s)
        needs_step && (fs[lf] = residual/s)
        needs_ramp && (fa[lf] = residual*(-expm1(-s*input.twindow))/s^2)
    end
    for lf in 2:nfg
        needs_step && (fs[ntg-lf+2] = conj(fs[lf]))
        needs_ramp && (fa[ntg-lf+2] = conj(fa[lf]))
    end
    factor = if fft_workspace === nothing
        tg = (0:ntg-1).*dtg
        dfg.*exp.(alpha.*tg)
    else
        fft_workspace.factor
    end
    step,ramp = Float64[],Float64[]
    if needs_step
        transformed = fft_plan === nothing ? bfft(fs) : fft_plan*fs
        step = fft_workspace === nothing ? real.(transformed).*factor :
            (fft_workspace.step .= real.(transformed).*factor)
    end
    if needs_ramp
        transformed = fft_plan === nothing ? bfft(fa) : fft_plan*fa
        ramp = fft_workspace === nothing ? real.(transformed).*factor :
            (fft_workspace.ramp .= real.(transformed).*factor)
    end
    # Set the step and ramp responses to zero at t=0 independently; this
    # avoids a -t*step(0) drift from integrating the truncated step response.
    needs_ramp && (ramp .-= ramp[1])
    needs_step && (step .-= step[1])
    function sample(values,lag)
        lag <= 0 && return 0.0
        x = lag/dtg
        i = min(floor(Int,x)+1,ntg)
        b = x-floor(x)
        (1-b)*values[i]+b*values[min(i+1,ntg)]
    end
    out = Vector{Float64}(undef,input.nt)
    dt = input.twindow/(input.nt-1)
    for it in eachindex(out)
        t = (it-1)*dt
        value = iszero(input.sinj0[1]) ? 0.0 : input.sinj0[1]*sample(step,t)
        for j in 2:length(input.ts0)
            t0,t1 = input.ts0[j-1],input.ts0[j]
            dh = input.sinj0[j]-input.sinj0[j-1]
            iszero(dh) && continue
            value += t1 == t0 ? dh*sample(step,t-t1) :
                dh/(t1-t0)*(sample(ramp,t-t0)-sample(ramp,t-t1))
        end
        iszero(linear) || (value += linear*source_slope(input,t))
        iszero(half) || (value += half*source_half_derivative(input,t))
        iszero(direct) || (value += direct*source_value(input,t))
        iszero(inverse_half) || (value += inverse_half*source_half_integral(input,t))
        out[it] = value
    end
    out
end
