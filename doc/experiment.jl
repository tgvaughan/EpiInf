module Phylodyn

function simSIS(S0, I0, beta, gamma, T)

    t = [0.0]
    S = [S0]
    I = [I0]

    while true
        aInfect = beta*S[end]*I[end]
        aRemove = gamma*I[end]
        a0 = aInfect + aRemove

        if a0>0
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

function drawCoalCount(Itraj, k)

    dI = diff(Itraj)
    I = Itraj[end]

    r = 0
    for i in length(dI):-1:1
        if dI[i]==1
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

    return k
end

end
