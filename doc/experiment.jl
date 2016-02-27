module Phylodyn

using Distributions

function simSIS(S0::Int, I0::Int, beta, gamma, T)

    t = [0.0]
    S = [S0]
    I = [I0]

    while true
        aInfect = beta*S[end]*I[end]
        aRemove = gamma*I[end]
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

    return t, I, S
end

function simB(I0, λ, T)
    t = [0.0]
    I = [I0]

    while true
        a0 = λ*I[end]

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

function mainB(;N=10, λ=1, k=10, T=1)
    iter = 10000

    r = zeros(Int, iter)
    rp = zeros(Int, iter)
    rpp = zeros(Int, iter)
    n = zeros(Int, iter)
    npp = zeros(Int, iter)
    for i in 1:iter

        t = []
        I = []
        while true
            t, I = simB(N, λ, T)

#            if I[end]>=k
                break
#            end
        end

        n[i] = sum(map(x->x==1, diff(I)))
        r[i] = drawCoalCount(I, k)

        M = N + n[i]
        p = binomial(k,2)/binomial(M,2)
        if p>0
            if p<1
                rp[i] = rand(Binomial(n[i], p))
            else
                rp[i] = n[i]
            end
        end

        p = binomial(k,2)/binomial(N+1,2)
        x = p*λ*N*T
        y = (1-p)*λ*N*T
        rpp[i] = rand(Poisson(x))
        npp[i] = rand(Poisson(y)) + rpp[i]
    end

    return n, r, rp, rpp, npp
end


function mainSIS(;S0=150, I0=50, β=0.05, γ=0.1, k=20, T=0.2)
    iter = 2000

    r = zeros(Int, iter)
    rp = zeros(Int, iter)
    n = zeros(Int, iter)
    np = zeros(Int, iter)
    for i in 1:iter

        t = []
        I = []
        while true
            t, I, S = simSIS(S0, I0, β, γ, T)

#            if I[end]>=k
                break
#            end
        end

        n[i] = sum(map(x->x==1, diff(I)))
        r[i] = drawCoalCount(I, k)

        M = I0 + n[i]
        p = binomial(k,2)/binomial(M,2)
        if (p>0) 
            if (p<1)
                rp[i] = rand(Binomial(n[i], p))
            else
                rp[i] = n[i]
            end
        end

        p = binomial(k,2)/binomial(I0+1,2)
        x = p*β*S0*I0*T
        y = (1-p)*β*S0*I0*T
        z = γ*I0*T
        rp[i] = rand(Poisson(x))
        np[i] = rand(Poisson(y)) + rp[i]
    end

    return n, r, np, rp
end


end
