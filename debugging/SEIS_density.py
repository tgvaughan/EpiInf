#!/usr/bin/env python

from argparse import ArgumentParser, FileType
from sys import argv, exit
import json
import scipy as sp

from Tree import *


class Params:
    def __init__(self, beta, alpha, gamma, psi, N):
        self.beta = beta
        self.alpha = alpha
        self.gamma = gamma
        self.psi = psi
        self.N = N

class ParticleState:
    def __init__(self, S, E, I, lineages):
        self.S, self.E, self.I = S, E, I
        self.lineages = lineages

    def updatePropensities(self, params):
        self.infectionProp = params.beta*self.I*self.S
        self.activationProp = params.alpha*self.E
        self.recoveryProp = params.gamma*self.I
        self.samplingProp = params.psi*self.I
        self.totalNonSamplingProp = self.infectionProp + self.activationProp + self.recoveryProp

    def replaceWith(self, other):
        self.S, self.E, self.I = other.S, other.E, other.I
        self.lineages = other.lineages

    def initTraj(self):
        return {'t': [], 'S': [], 'I': [], 'E': []}

    def updateTraj(self, t, traj):
        traj['t'].append(t)
        traj['S'].append(self.S)
        traj['E'].append(self.E)
        traj['I'].append(self.I)


def updateParticle(particleState, params, t0, node):
    """Updates a single particle by simulating its trajectory over
    a tree interval and computing its weight."""

    P = 1.0
    t = t0

    thisTraj = particleState.initTraj()

    while True:
        particleState.updateTraj(t, thisTraj)

        particleState.updatePropensities(params)

        if particleState.totalNonSamplingProp > 0.0:
            deltat = sp.random.exponential(scale=1.0/particleState.totalNonSamplingProp)
        else:
            deltat = float("inf")

        if t + deltat>node.time:
            deltat = node.time - t
            endReached = True
        else:
            endReached = False

        # Incorporate probability of no sampling into weight.
        P *= sp.exp(-deltat*particleState.samplingProp)

        t += deltat

        if endReached:
            break;

        u = sp.random.uniform(low=0.0, high=particleState.totalNonSamplingProp)

        # Infection
        u -= particleState.infectionProp
        if u < 0.0:
            particleState.S -= 1
            particleState.E += 1
            continue

        # Activation
        u -= particleState.activationProp
        if u < 0.0:
            particleState.E -= 1
            particleState.I += 1
            continue

        # Recovery
        u -= particleState.recoveryProp
        if u < 0.0:
            particleState.I -= 1
            particleState.S += 1
            continue

        raise Exception("Event selection fell through.")

    if node.isLeaf():
        # Sample
        P *= params.psi
    else:
        # Coalescence
        P *= 1
            
    return (P, thisTraj)


def computeLikelihood(tree, params, Nparticles=1000):
    """Computes the likelihood of the model parameters given the tree."""

    logP = 0.0

    # Get sorted list of tree nodes:
    def cmpFunc(n1, n2):
        if n1.time < n2.time:
            return -1
        if n1.time > n2.time:
            return 1
        return 0
    nodeList = sorted(tree.getNodes(), cmp=cmpFunc)

    # Initialize particles
    particles = []
    particlesPrime = []
    for pidx in range(Nparticles):
        particles.append(ParticleState(params.N-1, 0, 1,  {'I': [nodeList[0]]}))
        particlesPrime.append(ParticleState(params.N-1, 0, 1, {'I': [nodeList[0]]}))

    # Initialize weights
    weights = sp.ones(Nparticles)

    forJson = {'intervals': []}
    lastTime = 0.0
    for i in range(len(nodeList)):

        print nodeList[i].time
        forJson['intervals'].append([])

        # Update particles
        for pidx, pState in enumerate(particles):
            (weights[pidx], traj) = updateParticle(pState, params, lastTime, nodeList[i])
            forJson['intervals'][i].append(traj)

        # Update likelihood
        logP += sp.log(sp.mean(weights))

        weights = weights/sum(weights)

        # Resample particles
        for newParticle in particlesPrime:
            newParticle.replaceWith(sp.random.choice(particles, p=weights))
        particles, particlesPrime = particlesPrime, particles

        lastTime = nodeList[i].time

    return (logP, forJson)


### MAIN ###
if __name__ == '__main__':

    parser = ArgumentParser(description="Compute parameter likelihood of SEIS transmission tree.")
    parser.add_argument("treefile", type=FileType('r'), help="File containing transmission tree (nexus/newick format)")
    parser.add_argument("-b", type=float, dest='beta', default=0.01, help="Infection rate")
    parser.add_argument("-a", type=float, dest='alpha', default=0.1, help="Activation rate")
    parser.add_argument("-g", type=float, dest='gamma', default=0.1, help="Recovery rate")
    parser.add_argument("-s", type=float, dest='psi', default=0.1, help="Psi sampling rate")
    parser.add_argument("-N", type=int, dest='N', default=100, help="Population size")
    parser.add_argument("--particles", type=int, default=1000, help="Number of particles to use in PF algorithm.")
    parser.add_argument("--jsonOutputFile", type=FileType('w'), help="File to which particle trajectories are written"
            "in JSON format")

    args = parser.parse_args()

    params = Params(args.beta, args.alpha, args.gamma, args.psi, args.N)

    logP, forJson = computeLikelihood(Tree(args.treefile), params, Nparticles=args.particles)

    print logP

    if args.jsonOutputFile != None:
        json.dump(forJson, args.jsonOutputFile)

