/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package epiinf.distribs;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import beast.util.TreeParser;
import epiinf.*;
import epiinf.models.EpidemicModel;
import epiinf.models.SISModel;
import epiinf.util.ReplacementSampler;
import feast.fileio.TreeFromNewickFile;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static epiinf.util.EpiInfUtilityMethods.getLogOrientedPoissonDensity;
import static epiinf.util.EpiInfUtilityMethods.getLogPoissonProb;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
    + "parameters.  This variant leaps over tree events.")
@Citation("Gabriel Leventhal, Timothy Vaughan, David Welch, Alexei Drummond, Tanja Stadler,\n" +
        "\"Exact phylodynamic inference using particle filtering\", in preparation.")
public class LeapingSMCTreeDensity extends TreeDistribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    public Input<Double> epsilonInput = new Input<>(
            "tauLeapingEpsilon", "Relative fraction of propensity change to allow " +
            "when selecting leap size.", 0.03);

    public Input<Integer> nResamplesInput = new Input<>(
            "nResamples",
            "Number of evenly-spaced particle resampling times.",
            100);

    EpidemicModel model;
    TreeEventList treeEventList;
    int nParticles, nResamples;
    double epsilon;

    // Keep these around so we don't have to create these arrays/lists
    // for every density evaluation.

    double[] particleWeights;
    EpidemicState[] particleStates, particleStatesNew;

    double recordedOrigin;
    List<EpidemicState> recordedTrajectory;
    List<List<EpidemicState>> particleTrajectories, particleTrajectoriesNew;


    public LeapingSMCTreeDensity() {
        treeIntervalsInput.setRule(Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() throws Exception {
        model = modelInput.get();
        if (model.treeOriginInput.get() == null)
            throw new IllegalArgumentException("The treeOrigin input to " +
                    "EpidemicModel must be set when the model is used for inference.");
        treeEventList = new TreeEventList(treeInput.get(), model.treeOriginInput.get());
        nParticles = nParticlesInput.get();
        nResamples = nResamplesInput.get();
        epsilon = epsilonInput.get();

        particleWeights = new double[nParticles];
        particleStates = new EpidemicState[nParticles];
        particleStatesNew = new EpidemicState[nParticles];

        recordedTrajectory = new ArrayList<>();
        particleTrajectories = new ArrayList<>();
        particleTrajectoriesNew = new ArrayList<>();

        for (int p=0; p<nParticles; p++) {
            particleTrajectories.add(new ArrayList<>());
            particleTrajectoriesNew.add(new ArrayList<>());
        }
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = calculateLogP(false);
        return logP;
    }

    public double calculateLogP(boolean recordTrajectory) {

        double thisLogP = 0.0;

        if (recordTrajectory) {
            recordedOrigin = treeEventList.getOrigin();
            recordedTrajectory.clear();
        }

        // Early exit if first tree event occurs before origin.
        if (treeEventList.getEventList().get(0).time < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Initialize particles and trajectory storage
        for (int p = 0; p < nParticles; p++) {
            particleStates[p] = model.getInitialState();

            if (recordTrajectory) {
                particleTrajectories.get(p).clear();
                particleTrajectoriesNew.get(p).add(model.getInitialState());
            }
        }

        double dtResamp = model.getOrigin()/(nResamples-1);
        for (int ridx = 1; ridx<nResamples; ridx++) {
            double nextResampTime = ridx*dtResamp;

            // Update particles
            double sumOfWeights = 0.0;
            for (int p = 0; p < nParticles; p++) {

                double newLogWeight = recordTrajectory ?
                        updateParticle(particleStates[p], nextResampTime, particleTrajectories.get(p))
                        : updateParticle(particleStates[p], nextResampTime, null);

                double newWeight = Math.exp(newLogWeight);

                particleWeights[p] = newWeight;
                sumOfWeights += newWeight;
            }

            if (!(sumOfWeights > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }

            // Update marginal likelihood estimate
            thisLogP += Math.log(sumOfWeights / nParticles);

            // Normalize weights
            for (int i=0; i<nParticles; i++)
                particleWeights[i] = particleWeights[i]/sumOfWeights;

            // Sample particle with replacement
            ReplacementSampler replacementSampler = new ReplacementSampler(particleWeights);
            for (int p=0; p<nParticles; p++) {
                int srcIdx = replacementSampler.next();
                particleStatesNew[p] = particleStates[srcIdx].copy();

                if (recordTrajectory) {
                    particleTrajectoriesNew.get(p).clear();
                    particleTrajectoriesNew.get(p).addAll(particleTrajectories.get(srcIdx));
                }
            }

            // Switch particleStates and particleStatesNew
            EpidemicState[] tempStates = particleStates;
            particleStates = particleStatesNew;
            particleStatesNew = tempStates;

            // Switch particleTrajectories and particleTrajectoriesNew
            if (recordTrajectory) {
                List<List<EpidemicState>> tmpTrajs = particleTrajectories;
                particleTrajectories = particleTrajectoriesNew;
                particleTrajectoriesNew = tmpTrajs;
            }
        }

        // Choose arbitrary trajectory to log.
        if (recordTrajectory)
            recordedTrajectory.addAll(particleTrajectories.get(0));

        return thisLogP;
    }

    /**
     * Updates weight and state of particle.
     *
     * @param particleState State of particle
     * @param nextResampTime Time at which next bootstrap resampling will occur
     * @param particleTrajectory if non-null, add particle states to this trajectory
     *
     * @return log conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState, double nextResampTime,
                                  List<EpidemicState> particleTrajectory) {
        double conditionalLogP = 0;

        while (true) {

            model.calculatePropensities(particleState);

            // Determine length of proposed leap
            double infectionProp = model.propensities[EpidemicEvent.INFECTION];
            double removeProp = model.propensities[EpidemicEvent.RECOVERY]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE];
//            double tau = model.getTau(epsilon, particleState, infectionProp, removeProp);
            double tau = Double.POSITIVE_INFINITY;

            double nextModelEventTime = model.getNextModelEventTime(particleState);
            double trueDt = Math.min(tau, Math.min(nextModelEventTime, nextResampTime) - particleState.time);

            // Count tree events contained within leap interval
            int nInfections = 0;
            int nPsiSampleRemoves = 0;
            int nPsiSampleNoRemoves = 0;
            TreeEvent nextTreeEvent = treeEventList.getEventList()
                    .get(particleState.treeIntervalIdx);
            while (nextTreeEvent.time < particleState.time + trueDt
                    && !model.timesEqual(nextTreeEvent.time, nextModelEventTime)) {

                switch(nextTreeEvent.type) {
                    case COALESCENCE:
                        nInfections += 1;
                        break;

                    case LEAF:
                        nPsiSampleRemoves += 1;
                        break;

                    case SAMPLED_ANCESTOR:
                        nPsiSampleNoRemoves += 1;
                        break;
                }
                particleState.treeIntervalIdx += 1;
                nextTreeEvent = treeEventList.getEventList()
                        .get(particleState.treeIntervalIdx);
            }

            int k = treeEventList.getLineageCounts().get(particleState.treeIntervalIdx);
            double observedInfectProp = infectionProp*k*(k-1)/particleState.I/(particleState.I + 1);
            double unobservedInfectProp = infectionProp - observedInfectProp;

            model.incrementState(particleState,
                    EpidemicEvent.MultipleInfections(
                    (int) Math.round(Randomizer.nextPoisson(trueDt*unobservedInfectProp))
                            + nInfections));

            model.incrementState(particleState,
                    EpidemicEvent.MultipleRecoveries(
                     (int)Math.round(Randomizer.nextPoisson(
                             trueDt*model.propensities[EpidemicEvent.RECOVERY]))));

            model.incrementState(particleState,
                    EpidemicEvent.MultiplePsiSampleRemove(nPsiSampleRemoves));

            conditionalLogP += getLogOrientedPoissonDensity(observedInfectProp*trueDt, nInfections, trueDt)
                    + getLogOrientedPoissonDensity(model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]*trueDt, nPsiSampleRemoves, trueDt)
                    + getLogOrientedPoissonDensity(model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]*trueDt, nPsiSampleNoRemoves, trueDt);

            if (conditionalLogP == Double.NEGATIVE_INFINITY
                    || !particleState.isValid() || particleState.I < k)
                return Double.NEGATIVE_INFINITY;

            if (nextModelEventTime <= nextResampTime && particleState.time + tau > nextModelEventTime) {
                // Handle model events

                ModelEvent nextModelEvent = model.getNextModelEvent(particleState);
                if (nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING) {
                    // Handle rho-sampling events

                    if (!model.timesEqual(nextTreeEvent.time, nextModelEventTime))
                        return Double.NEGATIVE_INFINITY;

                    int I = (int) Math.round(particleState.I);
                    int nSamples = nextTreeEvent.multiplicity;
                    conditionalLogP += Binomial.logChoose(I, nSamples)
                            + nSamples*Math.log(nextModelEvent.rho)
                            + (I-nSamples)*Math.log(1.0 - nextModelEvent.rho);

                    particleState.treeIntervalIdx += 1;
                }

                particleState.time = nextModelEventTime;
                particleState.modelIntervalIdx += 1;

                if (model.timesEqual(nextModelEventTime, nextResampTime)) {
                    break;
                }

            } else {
                // Check for end of particle resampling interval

                if (particleState.time + tau > nextResampTime) {
                    particleState.time = nextResampTime;
                    break;
                } else {
                    particleState.time += tau;
                }
            }

            if (particleTrajectory != null)
                particleTrajectory.add(particleState.copy());

        }

        if (particleTrajectory != null)
            particleTrajectory.add(particleState.copy());

        if (!particleState.isValid())
            return Double.NEGATIVE_INFINITY; // Can occur due to susceptible pool depletion
        else
            return conditionalLogP;
    }

    public List<EpidemicState> getRecordedTrajectory() {
        return recordedTrajectory;
    }

    public double getRecordedOrigin() {
        return recordedOrigin;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    protected boolean requiresRecalculation() {
        treeEventList.makeDirty();
        return true;
    }

    @Override
    public void restore() {
        treeEventList.makeDirty();
        super.restore();
    }

    @Override
    public boolean isStochastic() {
        return true;
    }

    public static void main(String[] args) throws Exception {

//        TreeFromNewickFile tree = new TreeFromNewickFile();
//        tree.initByName("fileName", "examples/SIS_DensityMapLeap.tree.newick",
//                "IsLabelledNewick", true,
//                "adjustTipHeights", false);
        TreeParser tree = new TreeParser("(A:2,B:2):2;", false, false, true, 1);

        SISModel model = new SISModel();
        model.initByName("treeOrigin", new RealParameter("4.0"),
                "S0", new IntegerParameter("99"),
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingVariable", new RealParameter("0.0"),
                "removalProb", new RealParameter("1.0"),
                "rhoSamplingProb", new RealParameter("0.1"),
                "rhoSamplingTime", new RealParameter("4.0"));

        LeapingSMCTreeDensity density = new LeapingSMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "nParticles", 1000,
                "nResamples", 1001);

        SMCTreeDensity densityTrue = new SMCTreeDensity();
        densityTrue.initByName(
                "tree", tree,
                "model", model,
                "nParticles", 1000);

        System.out.println(density.calculateLogP(true));
        System.out.println(densityTrue.calculateLogP(true));

        System.out.println("t I S");
        for (EpidemicState state : density.recordedTrajectory) {
            System.out.println(state.time + " " + state.I + " " + state.S);
        }
    }

}
