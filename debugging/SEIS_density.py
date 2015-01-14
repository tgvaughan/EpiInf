#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import argv, exit

class TreeEvent:
    def __init__(self, age, isLeaf):
        self.age = age
        self.isLeaf = isLeaf
        self.time = None

    def __repr__(self):
        if self.isLeaf:
            nodeType = "Sample"
        else:
            nodeType = "Coalescence"

        return nodeType + " at age={}, t={}".format(self.age, self.time) 


def loadTreeEvents(treefile):
    treeEvents = []

    for line in treefile.readlines():
        lineElements = line.strip().split(" ")
        age = float(lineElements[0])
        isLeaf = (lineElements[1] == "0")
        treeEvents.append(TreeEvent(age, isLeaf))

    # Sort events in order of decreasing age
    treeEvents.sort(key=lambda event: event.age, reverse=True)

    # Add forward times
    for event in treeEvents:
        event.time = treeEvents[0].age - event.age

    return treeEvents


class Params:
    def __init__(self, beta, alpha, gamma, psi, N):
        self.beta = beta
        self.alpha = alpha
        self.gamma = gamma
        self.psi = psi
        self.N = N


class ParticleState:
    def __init__(self, S, E, I, R, kE, kI):
        self.S, self.E, self.I, self.R = S, E, I, R
        self.kE, self.KI = kE, kI

    def updatePropensities(self, params):
        self.infectionProp = params.beta*self.I*self.S
        self.activationProp = params.alpha*self.E
        self.recoveryProp = params.gamma*self.I
        self.samplingProp = params.psi*self.I
        self.totalNonSamplingProp = self.infectionProp + self.activationProp + self.recoveryProp

    def replaceWith(self, other):
        self.S, self.E, self.I, self.R = other.S, other.E, other.I, other.R
        self.kE, self.kI = other.kE, other.kI


def updateParticle(particleState, params, t0, finalTreeEvent):
    """Updates a single particle by simulating its trajectory over
    a tree interval and computing its weight."""

    P = 1.0
    t = t0

    while True:
        particleState.updatePropensities()

        if particleState.totalProp > 0.0:
            t += scipy.random.exponential(scale=1.0/particleState.totalProp)
        else:
            t += float("inf")

        if t>finalTreeEvent.time:
            break

        u = scipy.random.uniform(low=0.0, high=particleState.totalProp)

        # Infection
        u -= particleState.infectionProp
        if u < 0.0:
            continue

        # Activation
        u -= particleState.activationProp
        if u < 0.0:
            continue

        # Recovery
        u -= particleState.recoveryProp
        if u < 0.0:
            continue

    return P


def computeLikelihood(treeEvents, params, Nparticles=1000):
    """Computes the likelihood of the model parameters given the tree."""

    logP = 0.0

    # Initialize particles
    particles = []
    particlesPrime = []
    for pidx in range(Nparticles):
        particles.append(ParticleState(params.N-1, 0, 1, 0, 0, 1))
        particlesPrime.append(ParticleState(params.N-1, 0, 1, 0, 0, 1))

    # Initialize weights
    weights = ones(Nparticles)

    for i in range(1,len(treeEvents)):

        # Update particles
        for pidx, pState in enumerate(particles):
            weights[pidx] = updateParticle(pState, params, treeEvents[i-1].time, treeEvents[i])

        # Update likelihood
        logP += np.log(np.mean(weights))

        weights = weights/sum(weights)

        # Resample particles
        for newParticle in enumerate(particlesPrime):
            newParticle.replaceWith(scipy.random.choice(particles, p=weights))
        particles, particlesPrime = particlesPrime, particles

    return logP


### MAIN ###
if __name__ == '__main__':

    parser = ArgumentParser(description="Compute parameter likelihood of SEIS transmission tree.")
    parser.add_argument("treefile", type=FileType('r'), help="File containing transmission tree (ExpoTree format)")
    parser.add_argument("-b", type=float, dest='beta', default=0.01, help="Infection rate")
    parser.add_argument("-a", type=float, dest='alpha', default=0.1, help="Activation rate")
    parser.add_argument("-g", type=float, dest='gamma', default=0.1, help="Recovery rate")
    parser.add_argument("-s", type=float, dest='psi', default=0.1, help="Psi sampling rate")
    parser.add_argument("-N", type=int, dest='N', default=100, help="Population size")
    parser.add_argument("--particles", type=int, default=1000, help="Number of particles to use in PF algorithm.")

    args = parser.parse_args()

    eventList = loadTreeEvents(args.treefile)

    params = Params(args.beta, args.alpha, args.gamma, args.psi, args.N)

    print computeLikelihood(eventList, params)


