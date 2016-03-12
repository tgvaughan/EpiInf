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

import beast.core.*;
import beast.core.Input.Validate;
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import epiinf.*;
import epiinf.models.EpidemicModel;
import epiinf.util.IncidenceData;
import epiinf.util.ReplacementSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
    + "parameters.")
@Citation("Gabriel Leventhal, Timothy Vaughan, David Welch, Alexei Drummond, Tanja Stadler,\n" +
        "\"Exact phylodynamic inference using particle filtering\", in preparation.")
public class SMCTreeDensity extends TreeDistribution implements TrajectoryRecorder {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    public Input<Boolean> useTauLeapingInput = new Input<>(
            "useTauLeaping", "Whether to use tau leaping approximation.",
            false);

    public Input<Double> epsilonInput = new Input<>(
            "tauLeapingEpsilon", "Relative fraction of propensity change to allow " +
            "when selecting leap size.", 0.03);

    public Input<Double> resampThreshInput = new Input<>(
            "resampThresh",
            "Resampling performed when the effective relative number of " +
                    "particles drops below this threshold.",
            0.3);

    public Input<IncidenceData> incidenceDataInput = new Input<>(
            "incidenceData",
            "Times and numbers of unsequenced samples.");

    EpidemicModel model;
    ObservedEventsList observedEventsList;
    int nParticles;
    boolean useTauLeaping;
    double epsilon, resampThresh;

    // Keep these around so we don't have to create these arrays/lists
    // for every density evaluation.

    double[] logParticleWeights, particleWeights;
    EpidemicState[] particleStates, particleStatesNew;

    List<EpidemicState> recordedTrajectoryStates;
    List<List<EpidemicState>> particleTrajectories, particleTrajectoriesNew;


    public SMCTreeDensity() {
        treeIntervalsInput.setRule(Validate.FORBIDDEN);
    }

    @Override
    public void initAndValidate() {
        model = modelInput.get();
        if (model.treeOriginInput.get() == null)
            throw new IllegalArgumentException("The treeOrigin input to " +
                    "EpidemicModel must be set when the model is used for inference.");
        observedEventsList = new ObservedEventsList(treeInput.get(), incidenceDataInput.get(),
                model.treeOriginInput.get());
        nParticles = nParticlesInput.get();
        useTauLeaping = useTauLeapingInput.get();
        epsilon = epsilonInput.get();
        resampThresh = resampThreshInput.get();

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

    @Override
    public double calculateLogP() {
        logP = calculateLogP(false);
        return logP;
    }

    public double calculateLogP(boolean recordTrajectory) {

        double thisLogP = 0.0;

        if (recordTrajectory) {
            recordedTrajectoryStates.clear();
        }

        // Early exit if first tree event occurs before origin.
        if (observedEventsList.getEventList().get(0).time < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // Initialize particles and trajectory storage
        for (int p = 0; p < nParticles; p++) {
            particleStates[p].assignFrom(model.getInitialState());

            if (recordTrajectory) {
                particleTrajectories.get(p).clear();
                particleTrajectoriesNew.get(p).add(model.getInitialState());
            }
        }

        for (ObservedEvent observedEvent : observedEventsList.getEventList()) {

            // Update particles and record max log weight
            double maxLogWeight = Double.NEGATIVE_INFINITY;
            for (int p = 0; p < nParticles; p++) {

                logParticleWeights[p] += recordTrajectory ?
                        updateParticle(particleStates[p], particleTrajectories.get(p))
                        : updateParticle(particleStates[p], null);

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
                return Double.NEGATIVE_INFINITY;
            }

            double Neff = sumOfScaledWeights*sumOfScaledWeights/sumOfSquaredScaledWeights;

            if (Neff < resampThresh*nParticles || observedEvent.isFinal) {

                // Update marginal likelihood estimate
                thisLogP += Math.log(sumOfScaledWeights / nParticles) + maxLogWeight;

                // Normalize weights
                for (int i = 0; i < nParticles; i++)
                    particleWeights[i] = particleWeights[i] / sumOfScaledWeights;

                // Sample particle with replacement
                ReplacementSampler replacementSampler = new ReplacementSampler(particleWeights);
                for (int p = 0; p < nParticles; p++) {
                    int srcIdx = replacementSampler.next();
                    particleStatesNew[p].assignFrom(particleStates[srcIdx]);
                    logParticleWeights[p] = 0;

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
            recordedTrajectoryStates.addAll(particleTrajectories.get(0));

        return thisLogP;
    }

    /**
     * Updates weight and state of particle.
     *
     * @param particleState State of particle
     * @param particleTrajectory if non-null, add particle states to this trajectory
     *
     * @return log conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
                                  List<EpidemicState> particleTrajectory) {
        double conditionalLogP = 0;
        ObservedEvent nextObservedEvent = observedEventsList.getNextObservedEvent(particleState);

        while (true) {
            model.calculatePropensities(particleState);

            int lineages = nextObservedEvent.lineages;

            double infectionProp = model.propensities[EpidemicEvent.INFECTION];
            double unobservedInfectProp = infectionProp
                    *(1.0 - lineages * (lineages - 1) / particleState.I / (particleState.I + 1));
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

            boolean isLeap = particleTrajectory == null
                    && useTauLeaping;

            // Determine length of proposed leap and switch back to SSA
            // if length isn't much greater than the expected SSA step size.
            double tau = Double.POSITIVE_INFINITY;
            if (isLeap) {
                tau = model.getTau(epsilon, particleState, unobservedInfectProp, allowedRecovProp);
                if (tau < 10.0/allowedEventProp)
                    isLeap = false;
            }

            if (!isLeap) {

                // Determine size of time increment
                double dt;
                if (allowedEventProp > 0.0)
                    dt = Randomizer.nextExponential(allowedEventProp);
                else
                    dt = Double.POSITIVE_INFINITY;

                double nextModelEventTime = model.getNextModelEventTime(particleState);

                // Condition against psi-sampling and illegal recovery within interval
                double trueDt = Math.min(dt, Math.min(nextModelEventTime, nextObservedEvent.time) - particleState.time);
                conditionalLogP += -trueDt * (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                        + observedInfectProp + forbiddenRecovProp);

                // Increment time
                particleState.time += dt;

                // Deal with model events (rho sampling and rate shifts)
                if (nextModelEventTime < nextObservedEvent.time && particleState.time > nextModelEventTime) {

                    ModelEvent nextModelEvent = model.getNextModelEvent(particleState);
                    if (nextModelEvent.type == ModelEvent.Type.RHO_SAMPLING)
                        return Double.NEGATIVE_INFINITY;
                    else {
                        particleState.time = nextModelEventTime;
                        particleState.modelIntervalIdx += 1;
                        continue;
                    }
                }

                // Stop here if we're past the next observed event
                if (particleState.time > nextObservedEvent.time)
                    break;

                EpidemicEvent event = new EpidemicEvent();
                event.time = particleState.time;
                if (allowedEventProp * Randomizer.nextDouble() < unobservedInfectProp)
                    event.type = EpidemicEvent.INFECTION;
                else
                    event.type = EpidemicEvent.RECOVERY;


                model.incrementState(particleState, event);

                if (particleTrajectory != null)
                    particleTrajectory.add(particleState.copy());

                if (conditionalLogP == Double.NEGATIVE_INFINITY) {
                    // Should never get here, as we explicitly condition against
                    // events that cause this. However, rounding errors mock
                    // this rule.
                    return Double.NEGATIVE_INFINITY;
                }

            } else {

                double nextModelEventTime = model.getNextModelEventTime(particleState);
                double trueDt = Math.min(tau, Math.min(nextModelEventTime, nextObservedEvent.time) - particleState.time);
                conditionalLogP += -trueDt * (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                        + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                        + observedInfectProp + forbiddenRecovProp);

                EpidemicEvent infectEvent = new EpidemicEvent();
                infectEvent.type = EpidemicEvent.INFECTION;
                infectEvent.multiplicity = (int)Randomizer.nextPoisson(trueDt*unobservedInfectProp);

                EpidemicEvent recovEvent = new EpidemicEvent();
                recovEvent.type = EpidemicEvent.RECOVERY;
                recovEvent.multiplicity = (int)Randomizer.nextPoisson(trueDt*allowedRecovProp);

                model.incrementState(particleState, infectEvent);
                model.incrementState(particleState, recovEvent);

                if (conditionalLogP == Double.NEGATIVE_INFINITY
                        || !particleState.isValid() || particleState.I < lineages)
                    return Double.NEGATIVE_INFINITY;

                if (nextModelEventTime < nextObservedEvent.time && particleState.time + tau > nextModelEventTime) {
                    particleState.time = nextModelEventTime;
                    particleState.modelIntervalIdx += 1;
                } else {
                    if (particleState.time + tau > nextObservedEvent.time)
                        break;
                    else
                        particleState.time += tau;
                }
            }
        }

        particleState.time = nextObservedEvent.time;

        // Include probability of tree event
        if (nextObservedEvent.type == ObservedEvent.Type.COALESCENCE) {
            model.calculatePropensities(particleState);
            model.incrementState(particleState, EpidemicEvent.Infection);
            conditionalLogP += Math.log(2.0 / particleState.I / (particleState.I - 1)
                    * model.propensities[EpidemicEvent.INFECTION]);

        } else {

            if (model.timesEqual(nextObservedEvent.time, model.getNextModelEventTime(particleState))
                    && model.getNextModelEvent(particleState).type == ModelEvent.Type.RHO_SAMPLING) {

                ModelEvent nextModelEvent = model.getNextModelEvent(particleState);

                int I = (int) Math.round(particleState.I);
                int k = nextObservedEvent.multiplicity;
                conditionalLogP += Binomial.logChoose(I, k)
                        + k*Math.log(nextModelEvent.rho)
                        + (I-k)*Math.log(1.0 - nextModelEvent.rho);

                model.incrementState(particleState,
                        EpidemicEvent.MultipleRhoSamples(nextObservedEvent.multiplicity));

            } else {
                if (model.psiSamplingVariableInput.get() != null && nextObservedEvent.multiplicity == 1) {

                    model.calculatePropensities(particleState);

                    if (nextObservedEvent.type == ObservedEvent.Type.SAMPLED_ANCESTOR) {
                        conditionalLogP += Math.log(model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]/particleState.I);
                    } else {
                        double psiSamplingProp = (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                                + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]);

                        conditionalLogP += Math.log(psiSamplingProp);

                        boolean isRemoval = Randomizer.nextDouble()*psiSamplingProp
                                < model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE];

                        if (isRemoval) {
                            model.incrementState(particleState, EpidemicEvent.PsiSampleRemove);
                        } else {
                            conditionalLogP += Math.log(1.0 - (nextObservedEvent.lineages-1)/particleState.I);
                        }
                    }

                } else {
                    // No explicit sampling process
                    model.incrementState(particleState,
                            EpidemicEvent.MultipleOtherSamples(nextObservedEvent.multiplicity));
                }
            }

            conditionalLogP += GammaFunction.lnGamma(1 + nextObservedEvent.multiplicity);
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

        return new EpidemicTrajectory(null, recordedTrajectoryStates, observedEventsList.getOrigin());
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
        observedEventsList.makeDirty();
        return true;
    }

    @Override
    public void restore() {
        observedEventsList.makeDirty();
        super.restore();
    }

    @Override
    public boolean isStochastic() {
        return true;
    }

}
