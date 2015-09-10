package epiinf.distribs;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import beast.math.GammaFunction;
import beast.math.statistic.DiscreteStatistics;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Computes density of tree given stochastic SEIS model.")
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

    double[] aInfect, aActivate, aRecover;
    double aTotNS, aTotSamp;

    Particle[] particles, particlesPrime;
    double[] weights;

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

        aInfect = new double[tree.getNodeCount()];
        aActivate = new double[tree.getNodeCount()];
        aRecover = new double[tree.getNodeCount()];

        weights = new double[nParticles];
        particles = new Particle[nParticles];
        particlesPrime = new Particle[nParticles];
        for (int pIdx=0; pIdx<nParticles; pIdx++)
            particles[pIdx] = new Particle();
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        List<Node> treeNodes = new ArrayList<>(Arrays.asList(tree.getNodesAsArray()));
        treeNodes.sort((Node n1, Node n2) -> {
            if (n1.getHeight() > n2.getHeight())
                return -1;
            if (n1.getHeight() < n2.getHeight())
                return 1;
            return 0;
        });

        int S0 = S0Input.get().getValue();

        for (Particle particle : particles) {
            particle.S = S0;
            particle.E[0] = 1;
            particle.I[0] = 0;
        }

        for (int i=0; i<treeNodes.size(); i++) {

            double tStart = i>0
                    ? originInput.get().getValue() - treeNodes.get(i-1).getHeight()
                    : 0.0;

            // Skip zero-length intervals
//            if (originInput.get().getValue() - node.getHeight() == tStart)
//                continue;

            // Update particles
            for (int pIdx = 0; pIdx<nParticles; pIdx++)
                weights[pIdx] = updateParticle(particles[pIdx], treeNodes, i, tStart);

            // Compute average weight and include in tree density:
            logP += Math.log(DiscreteStatistics.mean(weights));

            // Particle resampling:
            for (int pIdx = 0; pIdx<particles.length; pIdx++) {
                int choice = Randomizer.randomChoicePDF(weights);
                particlesPrime[pIdx] =
                        particles[choice].getCopy();
            }

            Particle[] tmpParticles = particles;
            particles = particlesPrime;
            particlesPrime = tmpParticles;
        }

        // Account for arbitrary ordering of leaf nodes in treeNodes.
        logP += GammaFunction.lnGamma(nContemp+1);

        return logP;
    }

    protected double updateParticle(Particle particle,  List<Node> treeNodes, int intervalIdx, double tStart) {

        Node node = treeNodes.get(intervalIdx);

        double origin = originInput.get().getValue();
        double rateInfect = infectionRateInput.get().getValue();
        double rateActivate = activationRateInput.get().getValue();
        double rateRecover = recoveryRateInput.get().getValue();
        double rateSamp = psiSamplingRateInput.get().getValue();
        double pSamp = rhoSamplingProbInput.get().getValue();

        double t = tStart;
        double tEnd = origin - node.getHeight();

        double weight = 1.0;

        while(true) {

            // Calculate propensities

            aTotSamp = 0;
            aTotNS = 0;

            for (int j = 0; j <= intervalIdx; j++) {
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
                weight *= Math.exp(-(tEnd-t)*aTotSamp);
                break;
            }

            t += dt;
            weight *= Math.exp(-(tEnd-t)*aTotSamp);

            // Choose and implement reaction:
            double u = Randomizer.nextDouble()*aTotNS;

            for (int j=0; j<=intervalIdx; j++) {

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
        if (particle.I[intervalIdx] >= 1)
            weight *= particle.I[intervalIdx]/(double)(particle.E[intervalIdx] + particle.I[intervalIdx]);
        else
            return 0.0;

        // Compute probability of observed event:

        if (node.isLeaf()) {
            // Sampling event

            /* Following condition not exactly right as it is sometimes impossible
               to determine whether final sampling event was produced by rho
               sampling or by psi sampling. */
            if (rateSamp > 0.0 && (pSamp == 0.0 || node.getHeight()>0.0)) {
                weight *= aTotSamp;
            }

            if (pSamp > 0.0 && node.getHeight() == 0.0) {
                weight *= pSamp*Math.pow(1.0-pSamp, particle.I[intervalIdx]-1)
                        *particle.I[intervalIdx];
            }

            particle.I[intervalIdx] -= 1;
            particle.S += 1;

        } else {
            // Branching event

            int leftIdx = treeNodes.indexOf(node.getLeft());
            int rightIdx = treeNodes.indexOf(node.getRight());

            weight *= rateInfect;
            particle.S -= 1;
            if (Randomizer.nextBoolean()) {
                particle.E[leftIdx] = 1;
                particle.I[rightIdx] = 0;
            } else {
                particle.E[rightIdx] = 1;
                particle.I[leftIdx] = 0;
            }
        }

        return weight;
    }

}

