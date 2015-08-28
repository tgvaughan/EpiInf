package epiinf.distribs;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.math.Binomial;
import beast.math.statistic.DiscreteStatistics;
import beast.util.Randomizer;
import epiinf.EpidemicState;
import epiinf.models.SEISModel;

import java.util.*;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SEISTreeDensity extends TreeDistribution {

    public Input<IntegerParameter> S0Input = new Input<>(
            "S0", "Initial size of susceptible population.", Input.Validate.REQUIRED);

    public Input<RealParameter> infectionRateInput = new Input<>(
            "infectionRate", "Infection rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> activationRateInput = new Input<>(
            "activationRate", "Activation rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> recoveryRateInput = new Input<>(
            "recoveryRate", "Recovery rate.", Input.Validate.REQUIRED);

    public Input<RealParameter> psiSamplingRateInput = new Input<>(
            "psiSamplingRate",
            "Probability with which recoveries are translated into samples",
            Input.Validate.REQUIRED);

    public Input<RealParameter> rhoSamplingProbInput = new Input<>(
            "rhoSamplingProb",
            "Probability with which a lineage at the corresponding time"
                    + "is sampled.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin",
            "Origin of process.",
            Input.Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles",
            "Number of particles to use.",
            1000);

    TreeInterface tree;
    int nParticles;
    int nContemp;

    class Particle {
        int S;
        int[] E, I;

        public Particle() {
            S = 0;
            E = new int[tree.getNodeCount()];
            I = new int[tree.getNodeCount()];
        }

        public Particle getCopy() {
            Particle copy = new Particle();

            copy.S = S;
            System.arraycopy(E, 0, copy.E, 0, E.length);
            System.arraycopy(I, 0, copy.I, 0, I.length);

            return copy;
        }
    }

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        tree = treeInput.get();
        nParticles = nParticlesInput.get();

        // Count number of contemporaneously sampled leaves

        nContemp = 0;
        for (Node leaf : tree.getExternalNodes()) {
            if (leaf.getHeight()>0.0)
                continue;

            nContemp += 1;
        }
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        double origin = originInput.get().getValue();
        double rateInfect = infectionRateInput.get().getValue();
        double rateActivate = activationRateInput.get().getValue();
        double rateRecover = recoveryRateInput.get().getValue();
        double rateSamp = psiSamplingRateInput.get().getValue();
        double pSamp = rhoSamplingProbInput.get().getValue();

        List<Node> treeNodes = new ArrayList<>(Arrays.asList(tree.getNodesAsArray()));
        treeNodes.sort((Node n1, Node n2) -> {
            if (n1.getHeight() > n2.getHeight())
                return -1;
            if (n1.getHeight() < n2.getHeight())
                return 1;
            return 0;
        });

        int S0 = S0Input.get().getValue();

        Particle[] particles = new Particle[nParticles];
        Particle[] particlesPrime = new Particle[nParticles];
        double[] weights = new double[nParticles];

        for (int i=0; i<nParticles; i++) {
            Particle particle = new Particle();

            particle.S = S0;
            particle.E[0] = 1;
            particle.I[0] = 0;

            particles[i] = particle;
        }

        double[] aInfect = new double[treeNodes.size()];
        double[] aActivate = new double[treeNodes.size()];
        double[] aRecover = new double[treeNodes.size()];
        double aTotNS, aTotSamp;

        for (int i=0; i<treeNodes.size(); i++) {

            Node node = treeNodes.get(i);

            for (int pIdx = 0; pIdx<nParticles; pIdx++) {
                Particle particle = particles[pIdx];
                weights[i] = 1.0;

                double t;
                if (i > 0)
                    t = origin - treeNodes.get(i-1).getHeight();
                else
                    t = 0.0;

                double tEnd = origin - node.getHeight();

                while(true) {

                    // Calculate propensities

                    aTotSamp = 0;
                    aTotNS = 0;

                    for (int j = 0; j <= i; j++) {
                        aInfect[j] = rateInfect * particle.S * particle.I[j];
                        aTotNS += aInfect[j];

                        aActivate[j] = rateActivate * particle.E[j];
                        aTotNS += aActivate[j];

                        aRecover[j] = rateRecover * particle.I[j];
                        aTotNS += aRecover[j];

                        aTotSamp += rateSamp * particle.I[j];
                    }

//                    aTotNS = aInfectTot + aActivateTot + aRecoverTot;

                    // Increment time

                    double dt;
                    if (aTotNS>0.0)
                        dt = Randomizer.nextExponential(aTotNS);
                    else
                        dt = Double.POSITIVE_INFINITY;

                    if (t + dt>tEnd) {
                        weights[pIdx] *= Math.exp(-(tEnd-t)*aTotSamp);
                        break;
                    }

                    t += dt;
                    weights[pIdx] *= Math.exp(-(tEnd-t)*aTotSamp);

                    // Choose and implement reaction:
                    double u = Randomizer.nextDouble()*aTotNS;

                    for (int j=0; j<=i; j++) {

                        if (u < aInfect[j]) {
                            particle.S -= 1;
                            particle.E[j] += 1;
                            break;
                        }
                        u -= aInfect[j];

                        if (u < aActivate[j]) {
                            particle.E[j] -= 1;
                            particle.I[j] += 1;
                            break;
                        }
                        u -= aActivate[j];


                        if (u < aRecover[j]) {
                            particle.I[j] -= 1;
                            particle.S += 1;
                            break;
                        }
                        u -= aRecover[j];
                    }
                }

                // Compute probability of observed state (always I):
                weights[i] *= particle.I[i]/(double)(particle.E[i] + particle.I[i]);

                // Compute probability of observed event:

                if (node.isLeaf()) {
                    // Sampling event

                    if (rateSamp > 0.0) {
                        weights[i] *= rateSamp;
                        particle.I[i] -= 1;
                        particle.S += 1;
                    }

                } else {
                    // Branching event

                    weights[i] *= rateInfect;
                    particle.S -= 1;
                    particle.E[i] += 1;
                }
            }

            // Compute average weight and include in tree density:
            logP += Math.log(DiscreteStatistics.mean(weights));

            // Particle resampling:
            for (int pIdx = 0; pIdx<particles.length; pIdx++) {
                particlesPrime[pIdx] =
                        particles[Randomizer.randomChoicePDF(weights)].getCopy();
            }

            Particle[] tmpParticles = particles;
            particles = particlesPrime;
            particlesPrime = tmpParticles;
        }

        return logP;
    }


}

