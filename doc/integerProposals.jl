module IntProposalTests

using Distributions

function uniformTarget(n,x,y)
    if x<0.0 || x>1.0 || y<0.0 || y>1.0 || n<2 || n>100
        return -Inf
    else
        return 0.0
    end
end

function mcmc(;target=uniformTarget, iter=1000000)
    n = zeros(Int,iter)
    x = zeros(iter)
    y = zeros(iter)

    n[1] = rand(0:100)
    x[1] = rand(Uniform(0,1))
    y[1] = rand(Uniform(0,1))

    lastLogP = target(n[1], x[1], y[1])

    for i in 2:iter

        # Proposal
        prop = rand(1:4)

        logHR = 0.0

        if prop == 1
            np = n[i-1] + rand(-10:10)
            xp = x[i-1]
            yp = y[i-1]
        elseif prop == 2
            np = n[i-1]
            xp = x[i-1] + rand(Uniform(-0.1,0.1))
            yp = y[i-1]
        elseif prop == 3
            np = n[i-1]
            xp = x[i-1]
            yp = y[i-1] + rand(Uniform(-0.1,0.1))
        else
            #np = n[i-1] + rand(-10:10)
            np = rand(Geometric(1/n[i-1]))
            xp = x[i-1]*n[i-1]/np
            yp = y[i-1]*n[i-1]/np

            if np>1
                logHR = 2*log(n[i-1]/np) + log(pdf(Geometric(1/np),n[i-1])/pdf(Geometric(1/n[i-1]), np))
            end
        end

        alpha = target(np,xp,yp) - lastLogP + logHR

        if (alpha>0.0 || rand(Uniform(0,1))<exp(alpha)) && isfinite(xp) && isfinite(np)
            # accept

            n[i] = np
            x[i] = xp
            y[i] = yp
        else
            # reject

            n[i] = n[i-1]
            x[i] = x[i-1]
            y[i] = y[i-1]
        end
    end

    return n,x,y
end

end
