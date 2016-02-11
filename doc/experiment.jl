module Phylodyn

using Distributions

function simSIS(S0, I0, β, γ, T)

    t = [0.0]
    S = [S0]
    I = [I0]

    while true
        aInfect = β*S[end]*I[end]
        aRemove = γ*I[end]
        a0 = aInfect + aRemove

        if a0 > 0
            t = [t; t[end] + randexp()/a0]
        else
            t = [t; Inf]
        end

        if t[end] > T
            S = [S; S[end]]
            I = [I; I[end]]
            break
        end

        if rand()*a0 < aInfect
            # Infection
            S = [S; S[end]-1]
            I = [I; I[end]+1]
        else
            # Removal
            S = [S; S[end]+1]
            I = [I; I[end]-1]
        end
    end

    return t, S, I
end

function simB(I0, lambda, T)
    t = [0.0]
    I = [I0]

    while true
        a0 = lambda*I[end]

        if a0 > 0
            t = [t; t[end] + randexp()/a0]
        else
            t = [t; Inf]
        end

        if t[end] > T
            I = [I; I[end]]
            break
        end

        I = [I; I[end]+1]
    end

    return t, I
end

function drawCoalCount(Itraj, k)

    dI = reverse(diff(Itraj[1:end-1]))
    I = Itraj[end]

    r = 0
    for delta in dI
        if delta==1
            p = k*(k-1)/(I*(I-1))
            if rand()<p
                r += 1
                k -= 1
            end

            I -= 1
        else
            I += 1
        end
    end

    return r
end

function mainB(;N=10, k=10, T=1)
    iter = 10000

    r = zeros(iter)
    rp = zeros(iter)
    n = zeros(iter)
    for i in 1:iter

        t = []
        I = []
        while true
            t, I = simB(N, 1, T)

            if I[end]>=k
                break
            end
        end

        n[i] = sum(map(x->x==1, diff(I)))
        r[i] = drawCoalCount(I, k)

        M = N + n[i]
        p = k*(k-1)/M/(M-1)
        rp[i] = rand(Binomial(n[i], p))
    end

    return n, r, rp
end


function mainSIS(;S0=150, I0=50, β=0.05, γ=0.1, k=20, T=0.2)
    iter = 1000

    r = zeros(iter)
    rp = zeros(iter)
    rpp = zeros(iter)
    n = zeros(iter)
    for i in 1:iter

        t = []
        I = []
        while true
            t, I = simSIS(S0, I0, β, γ, T)

            if I[end]>=k
                break
            end
        end

        n[i] = sum(map(x->x==1, diff(I)))
        r[i] = drawCoalCount(I, k)

        M = I0 + n[i]
        p = k*(k-1)/M/(M-1)
        if (p>0 && p<1)
            rp[i] = rand(Binomial(n[i], p))
        end

        #p = k*(k-1)/(I0+1)/I0
        #rpp[i] = rand(Binomial(n[i], p))
    end

    return n, r, rp, rpp
end


end
