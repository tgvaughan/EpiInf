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

import beast.core.*;
import beast.core.Input.Validate;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import epiinf.*;
import epiinf.models.EpidemicModel;
import epiinf.util.EpiInfUtilityMethods;
import epiinf.util.ReplacementSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
    + "parameters.  This variant leaps over tree events.")
@Citation("Gabriel Leventhal, Timothy Vaughan, David Welch, Alexei Drummond, Tanja Stadler,\n" +
        "\"Exact phylodynamic inference using particle filtering\", in preparation.")
public class LeapingSMCTreeDensity extends EpiTreePrior {

    public Input<Function> incidenceInput = new Input<>(
            "incidenceParameter",
            "Ages of unsequenced samples.");

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    public Input<Double> epsilonInput = new Input<>(
            "tauLeapingEpsilon", "Relative fraction of propensity change to allow " +
            "when selecting leap size. A value of 1 causes leaps to bridge the " +
            "entire interval between resampling events.", 0.03);

    public Input<Integer> nResamplesInput = new Input<>(
            "nResamples",
            "Number of evenly-spaced particle resampling times. Whether or not " +
                    "resampling actually takes place at each of these times " +
                    "is determined by the particle weight variance via " +
                    "resampThresh.",
            100);

    public Input<Double> resampThreshInput = new Input<>(
            "resampThresh",
            "Resampling performed when the effective relative number of " +
                    "particles drops below this threshold.",
            0.5);

    int nParticles, nResamples;
    double epsilon, resampThresh;

    // Keep these around so we don't have to create these arrays/lists
    // for every density evaluation.

    double[] particleWeights, logParticleWeights;
    EpidemicState[] particleStates, particleStatesNew;

    List<EpidemicState> recordedEpidemicStates;
    List<List<EpidemicState>> particleTrajectories, particleTrajectoriesNew;


    public LeapingSMCTreeDensity() {
        treeIntervalsInput.setRule(Validate.FORBIDDEN);
        treeInput.setRule(Validate.OPTIONAL); // Possible to have only incidence data!
    }

    @Override
    public void initAndValidate() {
        model = modelInput.get();

        if (treeInput.get() == null && incidenceInput.get() == null)
            throw new IllegalArgumentException("Must specify at least one of tree or incidence.");

        observedEventsList = new ObservedEventsList(treeInput.get(), incidenceInput.get(),
                model.originInput.get(), finalSampleOffsetInput.get());
        nParticles = nParticlesInput.get();
        nResamples = nResamplesInput.get();
        epsilon = epsilonInput.get();
        resampThresh = resampThreshInput.get();

        particleWeights = new double[nParticles];
        logParticleWeights = new double[nParticles];
        particleStates = new EpidemicState[nParticles];
        particleStatesNew = new EpidemicState[nParticles];

        recordedEpidemicStates = new ArrayList<>();
        particleTrajectories = new ArrayList<>();
        particleTrajectoriesNew = new ArrayList<>();

        for (int p=0; p<nParticles; p++) {
            particleTrajectories.add(new ArrayList<>());
            particleTrajectoriesNew.add(new ArrayList<>());
        }
    }

    @Override
    public double calculateLogP() {
        logP = calculateLogP(false);
        return logP;
    }

    public double calculateLogP(boolean recordTrajectory) {

        double thisLogP = 0.0;

        if (recordTrajectory) {
            recordedEpidemicStates.clear();
        }

        // Early exit if first tree event occurs before origin.
        if (observedEventsList.getEventList().get(0).time < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Initialize particles and trajectory storage
        for (int p = 0; p < nParticles; p++) {
            particleStates[p] = model.getInitialState();
            logParticleWeights[p] = 0.0;

            if (recordTrajectory) {
                particleTrajectories.get(p).clear();
                particleTrajectoriesNew.get(p).add(model.getInitialState());
            }
        }

        double dtResamp = model.getOrigin()/(nResamples-1);
        for (int ridx = 1; ridx<nResamples; ridx++) {
            double nextResampTime = ridx*dtResamp;

            // Update particles
            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p = 0; p < nParticles; p++) {
                if (logParticleWeights[p] == Double.NEGATIVE_INFINITY)
                    continue;

                logParticleWeights[p] += recordTrajectory ?
                        updateParticle(particleStates[p], nextResampTime, particleTrajectories.get(p))
                        : updateParticle(particleStates[p], nextResampTime, null);

                maxLogWeight = Math.max(logParticleWeights[p], maxLogWeight);
            }

            // Compute mean of weights scaled relative to max log weight
            double sumOfScaledWeights = 0;
            double sumOfSquaredScaledWeights = 0;
            for (int p=0; p<nParticles; p++) {
                particleWeights[p] = Math.exp(logParticleWeights[p] - maxLogWeight);
                sumOfScaledWeights += particleWeights[p];
                sumOfSquaredScaledWeights += particleWeights[p]*particleWeights[p];
            }

            if (!(sumOfScaledWeights > 0.0)) {
//                System.out.println("Particle ensemble crashed.");
                return Double.NEGATIVE_INFINITY;
            }


            // Resample only when effective particle count is low

            double Neff = sumOfScaledWeights*sumOfScaledWeights/sumOfSquaredScaledWeights;
            if (Neff < resampThresh*nParticles || ridx == nResamples-1) {

                // Update marginal likelihood estimate
                thisLogP += Math.log(sumOfScaledWeights / nParticles) + maxLogWeight;
//                thisLogP += DiscreteStatistics.median(logParticleWeights);

                // Normalize weights
                for (int i = 0; i < nParticles; i++)
                    particleWeights[i] = particleWeights[i] / sumOfScaledWeights;

                // Sample particle with replacement
                ReplacementSampler replacementSampler = new ReplacementSampler(particleWeights);
                for (int p = 0; p < nParticles; p++) {
                    int srcIdx = replacementSampler.next();
                    particleStatesNew[p] = particleStates[srcIdx].copy();
                    logParticleWeights[p] = 0.0;

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
        }

        // Choose arbitrary trajectory to log.
        if (recordTrajectory)
            recordedEpidemicStates.addAll(particleTrajectories.get(0));

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
            double allowedRecovProp, forbiddenRecovProp;
            int k0 = observedEventsList.getCurrentLineageCount(particleState);
            if (particleState.I>k0) {
                allowedRecovProp = model.propensities[EpidemicEvent.RECOVERY];
                forbiddenRecovProp = 0.0;
            } else {
                allowedRecovProp = 0.0;
                forbiddenRecovProp = model.propensities[EpidemicEvent.RECOVERY];
            }
            double removeProp = allowedRecovProp
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE];

            double tau = epsilon < 1
                    ? model.getTau(epsilon, particleState, infectionProp, removeProp)
                    : Double.POSITIVE_INFINITY;

            double nextModelEventTime = model.getNextModelEventTime(particleState);
            double trueDt = Math.min(tau, Math.min(nextModelEventTime, nextResampTime) - particleState.time);

            // Waiting time weight contributions

            double unobservedInfectProp = infectionProp
                    *(1.0 - k0 * (k0 - 1) / particleState.I / (particleState.I + 1));
            double observedInfectProp = infectionProp - unobservedInfectProp;

            conditionalLogP += -trueDt*(model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                    + observedInfectProp + forbiddenRecovProp);

            // Count tree events contained within leap interval

            int nCoalescences = 0;
            int nLeaves = 0;
            int nSampledAncestors = 0;
            int nUnsequencedSamples = 0;
            ObservedEvent nextObservedEvent = observedEventsList.getEventList()
                    .get(particleState.observedEventIdx);
            while (model.timesLEQ(nextObservedEvent.time, particleState.time + trueDt)
                    && !model.timesEqual(nextObservedEvent.time, nextModelEventTime)) {

                switch(nextObservedEvent.type) {
                    case COALESCENCE:
                        nCoalescences += 1;
                        break;

                    case LEAF:
                        nLeaves += 1;
                        break;

                    case SAMPLED_ANCESTOR:
                        nSampledAncestors += 1;
                        break;

                    case UNSEQUENCED_SAMPLE:
                        nUnsequencedSamples += 1;
                        break;
                }
                particleState.observedEventIdx += 1;

                if (!nextObservedEvent.isFinal) {
                    nextObservedEvent = observedEventsList.getEventList()
                            .get(particleState.observedEventIdx);
                } else {
                    nextObservedEvent = null;
                    break;
                }
            }

            // Perform state increments

            int nUnobservedInfections = (int)Randomizer.nextPoisson(trueDt*unobservedInfectProp);
            int nInfections = nUnobservedInfections + nCoalescences;
            model.incrementState(particleState,
                    EpidemicEvent.MultipleInfections(nInfections));

            model.incrementState(particleState,
                    EpidemicEvent.MultipleRecoveries(
                            (int)Randomizer.nextPoisson(trueDt*allowedRecovProp)));

            int nSampleRemovals;
            if (nLeaves>0) {
                double sampleRemovalProb = model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] /
                        (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                                + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);
                nSampleRemovals = EpiInfUtilityMethods.nextBinomial(sampleRemovalProb, nLeaves);
                model.incrementState(particleState,
                        EpidemicEvent.MultiplePsiSampleRemove(nSampleRemovals));
            } else
                nSampleRemovals = 0;

            int nUnsequencedSampleRemovals;
            if (nUnsequencedSamples>0) {
                double sampleRemovalProb = model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] /
                        (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                                + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);
                nUnsequencedSampleRemovals = EpiInfUtilityMethods.nextBinomial(sampleRemovalProb, nUnsequencedSamples);
                model.incrementState(particleState,
                        EpidemicEvent.MultiplePsiSampleRemove(nUnsequencedSampleRemovals));

            } else {
                nUnsequencedSampleRemovals = 0;
            }

            // Compute weight contributions of tree events:

            int k = observedEventsList.getCurrentLineageCount(particleState);
            if (nCoalescences>0) {
                if (particleState.I<2)
                    return Double.NEGATIVE_INFINITY;

                conditionalLogP += nCoalescences*Math.log(
                        2.0/particleState.I/(particleState.I-1)
                        *model.propensities[EpidemicEvent.INFECTION]);
            }

            if (nSampledAncestors>0) {
                if (particleState.I<1)
                    return Double.NEGATIVE_INFINITY;

                conditionalLogP += nSampledAncestors*Math.log(
                        model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]/particleState.I);
            }

            if (nLeaves>0) {
                conditionalLogP += nLeaves*Math.log(
                        model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);

                if (nSampleRemovals<nLeaves) {
                    if (particleState.I<=k)
                        return Double.NEGATIVE_INFINITY;

                    conditionalLogP += (nLeaves - nSampleRemovals) * Math.log(1.0 - k / particleState.I);
                }
            }

            if (nUnsequencedSamples>0) {
                conditionalLogP += nUnsequencedSamples*Math.log(
                        model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);
            }

            // Discard if we're left with an invalid state

            if (conditionalLogP == Double.NEGATIVE_INFINITY
                    || !particleState.isValid() || particleState.I < k)
                return Double.NEGATIVE_INFINITY;

            // Do anything special that needs to be done at the end of the leap

            if (model.timesLEQ(nextModelEventTime, nextResampTime) && particleState.time + tau > nextModelEventTime) {

                // Handle model events

                ModelEvent nextModelEvent = model.getNextModelEvent(particleState);
                if (nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING) {
                    // Handle rho-sampling events

                    if (nextObservedEvent == null || !model.timesEqual(nextObservedEvent.time, nextModelEventTime))
                        return Double.NEGATIVE_INFINITY;

                    int I = (int) Math.round(particleState.I);
                    int nSamples = nextObservedEvent.multiplicity;
                    conditionalLogP += Binomial.logChoose(I, nSamples)
                            + nSamples*Math.log(nextModelEvent.rho)
                            + (I-nSamples)*Math.log(1.0 - nextModelEvent.rho)
                            + GammaFunction.lnGamma(1 + nSamples);

                    particleState.observedEventIdx += 1;
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

    @Override
    public EpidemicTrajectory getConditionedTrajectory() {
        calculateLogP(true);
        return new EpidemicTrajectory(null, recordedEpidemicStates, observedEventsList.getOrigin());
    }

}
