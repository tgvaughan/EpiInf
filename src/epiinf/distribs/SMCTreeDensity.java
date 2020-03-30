/*
* Copyright (C) 2013 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.math.statistic.DiscreteStatistics;
import beast.util.Randomizer;
import epiinf.*;
import epiinf.models.EpidemicModel;
import epiinf.util.ReplacementSampler;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
    + "parameters.")
@Citation("Timothy Vaughan, Gabriel Leventhal, David Welch, Alexei Drummond, Tanja Stadler,\n" +
        "\"Estimating epidemic incidence and prevalence curves from genomic data\",\n" +
        "Molecular Biology and Evolution, 36(8):1804-1816 (2019).")
public class SMCTreeDensity extends TreeDistribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Input.Validate.REQUIRED);

    public Input<Function> finalTreeSampleOffsetInput = new Input<>(
            "finalTreeSampleOffset",
            "Difference in time between final TREE sample and end " +
                    "of observation process.  (Defaults to zero.)");

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    public Input<Boolean> useTauLeapingInput = new Input<>(
            "useTauLeaping", "Whether to use tau leaping approximation.",
            false);

    public Input<Double> epsilonInput = new Input<>(
            "tauLeapingEpsilon", "Relative fraction of propensity change to allow " +
            "when selecting leap size.", 0.03);

    public Input<Integer> minLeapCountInput = new Input<>(
            "minLeapCount", "This is the minimum number of identically-sized " +
            "tau leaps that will be performed across the tree.", 100);

    public Input<Double> relStdThreshInput = new Input<>(
            "relStdThresh", "Threshold on relative size of standard deviation of" +
            "reaction firings below which deterministic approximation will be used. " +
            "Set to zero to turn off.", 0.0);

    public Input<Double> resampThreshInput = new Input<>(
            "resampThresh",
            "Resampling performed when the effective relative number of " +
                    "particles drops below this threshold.",
            0.3);

    public Input<RealParameter> incidenceParamInput = new Input<>(
            "incidenceParameter",
            "Ages of unsequenced samples.");

    public Input<IncidenceData> incidenceDataInput = new Input<>(
            "incidenceData",
            "Incindence data (as a histogram).");


    protected EpidemicModel model;
    protected ObservedEventsList observedEventsList;

    int nParticles;
    boolean useTauLeaping;
    double epsilon, resampThresh, relStdThresh;
    int minLeapCount;

    // Keep these around so we don't have to create these arrays/lists
    // for every density evaluation.

    double[] logParticleWeights, particleWeights;
    EpidemicState[] particleStates, particleStatesNew;

    List<EpidemicState> recordedTrajectoryStates;
    EpidemicTrajectory storedTrajectory;
    List<List<EpidemicState>> particleTrajectories, particleTrajectoriesNew;


    public SMCTreeDensity() {
        treeIntervalsInput.setRule(Validate.FORBIDDEN);
        treeInput.setRule(Validate.OPTIONAL); // Possible to have only incidence data!
    }

    @Override
    public void initAndValidate() {
        model = modelInput.get();

        if (treeInput.get() == null && incidenceParamInput.get() == null)
            throw new IllegalArgumentException("Must specify at least one of tree or incidence.");

        observedEventsList = new ObservedEventsList(treeInput.get(),
                incidenceParamInput.get(), incidenceDataInput.get(),
                model, finalTreeSampleOffsetInput.get());

        nParticles = nParticlesInput.get();

        useTauLeaping = useTauLeapingInput.get();
        epsilon = epsilonInput.get();
        minLeapCount = minLeapCountInput.get();
        resampThresh = resampThreshInput.get();
        relStdThresh = relStdThreshInput.get();

        particleWeights = new double[nParticles];
        logParticleWeights = new double[nParticles];
        particleStates = new EpidemicState[nParticles];
        particleStatesNew = new EpidemicState[nParticles];

        recordedTrajectoryStates = new ArrayList<>();
        particleTrajectories = new ArrayList<>();
        particleTrajectoriesNew = new ArrayList<>();

        for (int p=0; p<nParticles; p++) {
            particleTrajectories.add(new ArrayList<>());
            particleTrajectoriesNew.add(new ArrayList<>());

            particleStates[p] = new EpidemicState();
            particleStatesNew[p] = new EpidemicState();
        }
    }

    public double calculateLogP() {

        logP = 0.0;

        recordedTrajectoryStates.clear();

        // Early exit if first tree event occurs before origin.
        if (observedEventsList.getEventList().get(0).time < 0) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        // Initialize particles and trajectory storage
        for (int p = 0; p < nParticles; p++) {
            particleStates[p].assignFrom(model.getInitialState());
            logParticleWeights[p] = 0.0;

            particleTrajectoriesNew.get(p).clear();
            particleTrajectories.get(p).clear();
            particleTrajectories.get(p).add(model.getInitialState());
        }

        for (ObservedEvent observedEvent : observedEventsList.getEventList()) {
            if (!propagateEnsemble(observedEvent)) {
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
        }

        // Choose arbitrary trajectory to log.
        recordedTrajectoryStates.addAll(particleTrajectories.get(0));

        return logP;
    }

    /**
     * Propagate particle ensemble up to chosen observed event.
     *
     * @param nextObservedEvent Next observed event.
     *
     * @return true if propagation succeeds, false if it fails due to ensemble extinction
     */
    private boolean propagateEnsemble(ObservedEvent nextObservedEvent) {

            // Update particles and record max log weight
            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p = 0; p < nParticles; p++) {

                if (logParticleWeights[p] > Double.NEGATIVE_INFINITY)
                    logParticleWeights[p] += updateParticle(particleStates[p], particleTrajectories.get(p),
                            nextObservedEvent);

                maxLogWeight = Math.max(logParticleWeights[p], maxLogWeight);
            }


            // Compute mean of weights scaled relative to max log weight
            double sumOfScaledWeights = 0, sumOfSquaredScaledWeights = 0;
            for (int p=0; p<nParticles; p++) {
                particleWeights[p] = Math.exp(logParticleWeights[p] - maxLogWeight);
                sumOfScaledWeights += particleWeights[p];
                sumOfSquaredScaledWeights += particleWeights[p]*particleWeights[p];
            }

            if (!(sumOfScaledWeights > 0.0)) {
                return false;
            }

            double Neff = sumOfScaledWeights*sumOfScaledWeights/sumOfSquaredScaledWeights;

            if (Neff < resampThresh*nParticles || nextObservedEvent.type == ObservedEvent.Type.OBSERVATION_END) {
                // Update marginal likelihood estimate
                logP += Math.log(sumOfScaledWeights / nParticles) + maxLogWeight;

                // Normalize weights
                for (int i = 0; i < nParticles; i++)
                    particleWeights[i] = particleWeights[i] / sumOfScaledWeights;

                resampleParticles();
            }

            return true;
    }

    /**
     * Resample particle states from weighted particle distribution.
     */
    private void resampleParticles() {

        // Sample particle with replacement
        double[] tmpWeights = new double[nParticles];
        System.arraycopy(particleWeights, 0, tmpWeights, 0, nParticles);
        ReplacementSampler replacementSampler = new ReplacementSampler(tmpWeights);
        for (int p = 0; p < nParticles; p++) {
            int srcIdx = replacementSampler.next();
            particleStatesNew[p].assignFrom(particleStates[srcIdx]);
            logParticleWeights[p] = 0;

            particleTrajectoriesNew.get(p).clear();
            particleTrajectoriesNew.get(p).addAll(particleTrajectories.get(srcIdx));
        }

        // Switch particleStates and particleStatesNew
        EpidemicState[] tempStates = particleStates;
        particleStates = particleStatesNew;
        particleStatesNew = tempStates;

        // Switch particleTrajectories and particleTrajectoriesNew
        List<List<EpidemicState>> tmpTrajs = particleTrajectories;
        particleTrajectories = particleTrajectoriesNew;
        particleTrajectoriesNew = tmpTrajs;
    }

    /**
     * Updates weight and state of particle, simulating until the next tree
     * event (if there is one) or the last incidence report (if there isn't).
     *
     * @param particleState State of particle
     * @param particleTrajectory if non-null, add particle states to this trajectory
     * @param nextObservedEvent next observed event - termination of particle sim
     *
     * @return log conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
                                  List<EpidemicState> particleTrajectory,
                                  ObservedEvent nextObservedEvent) {
        double conditionalLogP = 0;
        ModelEvent nextModelEvent;
        double nextModelEventTime;
        double nextObservedEventTime = nextObservedEvent.time;

        double maxLeapSize = model.getOrigin()/minLeapCount;

        while (true) {
            nextModelEvent = model.getNextModelEvent(particleState);
            nextModelEventTime = model.getNextModelEventTime(particleState);

            model.calculatePropensities(particleState);

            int lineages = nextObservedEvent.lineages;

            double infectionProp = model.propensities[EpidemicEvent.INFECTION];
            double unobservedInfectProp = particleState.I > 0
                    ? infectionProp *(1.0 - lineages * (lineages - 1) / particleState.I / (particleState.I + 1))
                    : 0.0;
            double observedInfectProp = infectionProp - unobservedInfectProp;

            double allowedRecovProp, forbiddenRecovProp;
            if (particleState.I > lineages) {
                allowedRecovProp = model.propensities[EpidemicEvent.RECOVERY];
                forbiddenRecovProp = 0.0;
            } else {
                allowedRecovProp = 0.0;
                forbiddenRecovProp = model.propensities[EpidemicEvent.RECOVERY];
            }

            double allowedEventProp = unobservedInfectProp + allowedRecovProp;

            // Do we leap?

            boolean isLeap = useTauLeaping;

            // Determine length of proposed leap and switch back to SSA
            // if length isn't much greater than the expected SSA step size.
            double tau = maxLeapSize;
            if (isLeap) {
                if (epsilon>0.0) {
                    tau = Math.min(maxLeapSize,
                            model.getTau(epsilon, particleState, unobservedInfectProp, allowedRecovProp));

                }
            }

            if (tau < 10.0/allowedEventProp)
                isLeap = false;

            if (!isLeap) {
                particleState.algorithm = EpidemicState.Algorithm.SSA;

                // Determine size of time increment
                double dt;
                if (allowedEventProp > 0.0)
                    dt = Randomizer.nextExponential(allowedEventProp);
                else
                    dt = Double.POSITIVE_INFINITY;

                // Condition against psi-sampling and illegal recovery within interval
                double trueDt = Math.min(dt, Math.min(nextModelEventTime, nextObservedEventTime) - particleState.time);
                conditionalLogP += -trueDt * (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                        + observedInfectProp + forbiddenRecovProp);

                // Increment time
                particleState.time += dt;

                // Deal with model events (rho sampling and rate shifts)
                if (nextModelEventTime < nextObservedEventTime && particleState.time > nextModelEventTime) {

                    if (nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING) {
                        // Probability that rho sampling produced no samples
                        conditionalLogP += particleState.I*Math.log(1.0-nextModelEvent.rho);
                    }
                    particleState.time = nextModelEvent.time;
                    particleState.modelIntervalIdx += 1;
                    continue;
                }

                // Stop here if we're past the next observed event
                if (particleState.time > nextObservedEventTime)
                        break;

                EpidemicEvent event = new EpidemicEvent();
                event.time = particleState.time;
                if (allowedEventProp * Randomizer.nextDouble() < unobservedInfectProp)
                    event.type = EpidemicEvent.INFECTION;
                else
                    event.type = EpidemicEvent.RECOVERY;

                model.incrementState(particleState, event);

                if (conditionalLogP == Double.NEGATIVE_INFINITY) {
                    // Should never get here, as we explicitly condition against
                    // events that cause this. However, rounding errors mock
                    // this rule.
                    return Double.NEGATIVE_INFINITY;
                }

            } else {
                particleState.algorithm = EpidemicState.Algorithm.TL;

                double trueDt = Math.min(tau, Math.min(nextModelEventTime, nextObservedEventTime) - particleState.time);
                conditionalLogP += -trueDt * (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                        + observedInfectProp + forbiddenRecovProp);

                double propThresh = Double.POSITIVE_INFINITY;
                if (trueDt > 0.0 && relStdThresh > 0.0) {
                    if (relStdThresh < 1.0)
                        propThresh = 1.0/trueDt/relStdThresh/relStdThresh;
                    else
                        propThresh = 0.0;
                }

                EpidemicEvent infectEvent = new EpidemicEvent();
                infectEvent.type = EpidemicEvent.INFECTION;
                if (unobservedInfectProp<propThresh)
                    infectEvent.multiplicity = (int)Randomizer.nextPoisson(trueDt*unobservedInfectProp);
                else
                    infectEvent.multiplicity = (int)Math.round(trueDt*unobservedInfectProp);

                EpidemicEvent recovEvent = new EpidemicEvent();
                recovEvent.type = EpidemicEvent.RECOVERY;
                if (allowedEventProp < propThresh)
                    recovEvent.multiplicity = (int)Randomizer.nextPoisson(trueDt*allowedRecovProp);
                else
                    recovEvent.multiplicity = (int)Math.round(trueDt*allowedEventProp);

                model.incrementState(particleState, infectEvent);
                model.incrementState(particleState, recovEvent);

                if (conditionalLogP == Double.NEGATIVE_INFINITY
                        || !particleState.isValid() || particleState.I < lineages)
                    return Double.NEGATIVE_INFINITY;

                if (nextModelEventTime < nextObservedEvent.time && particleState.time + tau > nextModelEventTime) {
                    if (nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING) {
                        // Probability that rho sampling produces no samples
                        conditionalLogP += particleState.I*Math.log(1.0-nextModelEvent.rho);
                    }
                    particleState.time = nextModelEventTime;
                    particleState.modelIntervalIdx += 1;
                    continue;
                }

                // Stop here if we're at an observed event
                if (particleState.time + tau > nextObservedEventTime)
                    break;

                particleState.time += tau;
            }

            if (particleTrajectory != null)
                particleTrajectory.add(particleState.copy());
        }

        // Include probability of tree event and increment state if necessary
        if (nextObservedEvent.type != ObservedEvent.Type.OBSERVATION_END) {
            particleState.time = nextObservedEvent.time;
            conditionalLogP += getObservedEventProbability(particleState,
                    nextObservedEvent, nextObservedEventTime,
                    nextModelEvent, nextModelEventTime);
        } else
            particleState.time = model.getOrigin();

        if (particleTrajectory != null)
            particleTrajectory.add(particleState.copy());

        if (!particleState.isValid())
            return Double.NEGATIVE_INFINITY; // Can occur due to susceptible pool depletion
        else {
            particleState.observedEventIdx += 1;
            return conditionalLogP;
        }
    }

    private double getObservedEventProbability(EpidemicState particleState,
                                               ObservedEvent nextObservedEvent, double nextObservedEventTime,
                                               ModelEvent nextModelEvent, double nextModelEventTime) {

        double conditionalLogP = 0.0;

         if (nextObservedEvent.type == ObservedEvent.Type.COALESCENCE) {

             for (int i=0; i<nextObservedEvent.multiplicity; i++) {
                 model.calculatePropensities(particleState);

                 if (model.propensities[EpidemicEvent.INFECTION] == 0.0)
                     return Double.NEGATIVE_INFINITY;

                 model.incrementState(particleState, EpidemicEvent.Infection);
                 conditionalLogP += Math.log(2.0 / particleState.I / (particleState.I - 1)
                         * model.propensities[EpidemicEvent.INFECTION]);
             }

        } else {

            if (model.timesEqual(nextObservedEventTime, nextModelEventTime)
                    && nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING) {

                int I = (int) Math.round(particleState.I);
                int k = nextObservedEvent.multiplicity;
                conditionalLogP += Binomial.logChoose(I, k)
                        + k*Math.log(nextModelEvent.rho)
                        + (I-k)*Math.log(1.0 - nextModelEvent.rho);

                conditionalLogP += GammaFunction.lnGamma(1 + k);

                model.incrementState(particleState,
                        EpidemicEvent.MultipleRhoSamples(nextObservedEvent.multiplicity));

            } else {
                if (model.psiSamplingVariableInput.get() != null) {

                    // This awkward conditioning is to avoid unnecessary recalculation of
                    // reaction propensities.
                    if (nextObservedEvent.type == ObservedEvent.Type.UNSEQUENCED_SAMPLE)
                        model.calculatePropensities(particleState);

                    if (nextObservedEvent.type == ObservedEvent.Type.UNSEQUENCED_SAMPLE
                            && model.currentRemovalProb == 0.0) {

                        // Only works when removal disabled!!

                        if (particleState.I == 0) {
                            conditionalLogP = Double.NEGATIVE_INFINITY;
                        } else {
                            model.calculatePropensities(particleState);
                            conditionalLogP += nextObservedEvent.multiplicity * Math.log(model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);
                        }

                    } else {
                        for (int i = 0; i < nextObservedEvent.multiplicity; i++) {
                            model.calculatePropensities(particleState);

                            if (particleState.I == 0) {
                                conditionalLogP = Double.NEGATIVE_INFINITY;
                                break;
                            }

                            if (nextObservedEvent.type == ObservedEvent.Type.SAMPLED_ANCESTOR) {
                                conditionalLogP += Math.log(model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE] / particleState.I);
                            } else {
                                double psiSamplingProp = (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);

                                conditionalLogP += Math.log(psiSamplingProp);

                                boolean isRemoval;
                                if (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] == psiSamplingProp) {
                                    isRemoval = true;
                                } else {
                                    if (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE] == 0.0) {
                                        isRemoval = false;
                                    } else {
                                        isRemoval = Randomizer.nextDouble() * psiSamplingProp
                                                < model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE];
                                    }
                                }

                                if (isRemoval) {
                                    model.incrementState(particleState, EpidemicEvent.PsiSampleRemove);
                                } else {
                                    if (nextObservedEvent.type == ObservedEvent.Type.LEAF)
                                        conditionalLogP += Math.log(1.0 - (nextObservedEvent.lineages - 1) / particleState.I);
                                }
                            }

                            // Account for probability of sequencing (if non-null)
                            if (model.sequencingProbInput.get() != null) {
                                double seqProb = model.sequencingProbInput.get().getArrayValue();
                                if (nextObservedEvent.type == ObservedEvent.Type.UNSEQUENCED_SAMPLE)
                                    conditionalLogP += Math.log(1.0 - seqProb);
                                else
                                    conditionalLogP += Math.log(seqProb);
                            }
                        }
                    }

                } else {
                    // No explicit sampling process

                    // Note that we always need to condition on whether a sampling event was
                    // a removal or not.  Since this information is not provided in
                    // the data, we need to model at least this part of the sampling
                    // process.

                    for (int i=0; i<nextObservedEvent.multiplicity; i++) {
                        model.calculatePropensities(particleState);

                        if (nextObservedEvent.type == ObservedEvent.Type.SAMPLED_ANCESTOR) {
                            conditionalLogP += Math.log((1.0 - model.currentRemovalProb) / particleState.I);
                        } else {
                            boolean isRemoval = Randomizer.nextDouble() < model.currentRemovalProb;

                            if (isRemoval) {
                                model.incrementState(particleState, EpidemicEvent.OtherSampleRemove);
                            } else {
                                if (nextObservedEvent.type == ObservedEvent.Type.LEAF)
                                    conditionalLogP += Math.log(1.0 - (nextObservedEvent.lineages - 1)/particleState.I);
                            }
                        }

                        // Must explicitly check for validity of post-sampled state
                        // because there is no sampling rate/probability to prevent
                        // sampling when there are no lineages to sample.
                        if (!particleState.isValid())
                            return Double.NEGATIVE_INFINITY;
                    }
                }
            }

        }


        if (Double.isNaN(conditionalLogP)) {
            System.out.println("Particle weight is NaN - should not happen!");
        }

        return conditionalLogP;
    }

    /**
     * @return Trajectory conditioned on tree, for logging during PMMH analysis.
     */
    public EpidemicTrajectory getConditionedTrajectory() {
        return storedTrajectory;
    }

    /**
     * @return Most recent trajectory, for logging during trajectory mapping.
     */
    public EpidemicTrajectory getMostRecentTrajectory() {
        model.makeDirty();
        observedEventsList.makeDirty();
        calculateLogP();
        return new EpidemicTrajectory(null,
                new ArrayList<>(recordedTrajectoryStates),
                observedEventsList.getOrigin());
    }

    /**
     * @return Epidemic model
     */
    public EpidemicModel getModel() {
        return modelInput.get();
    }


    /*
     * StateNode implementation
     */

    @Override
    protected boolean requiresRecalculation() {
        observedEventsList.makeDirty();
        return true;
    }

    @Override
    public void restore() {
        observedEventsList.makeDirty();
        super.restore();
    }

    @Override
    protected void accept() {
        super.accept();

        List<EpidemicState> stateListCopy = new ArrayList<>(recordedTrajectoryStates);
        storedTrajectory = new EpidemicTrajectory(null, stateListCopy, observedEventsList.getOrigin());
    }

    /*
     * Distribution implementation
     */

    @Override
    public boolean isStochastic() {
        return true;
    }
}
