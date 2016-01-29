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
import beast.evolution.tree.TreeDistribution;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import epiinf.*;
import epiinf.models.EpidemicModel;
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

            // Update particles
            double sumOfWeights = 0.0;
            for (int p = 0; p < nParticles; p++) {

                double newLogWeight = recordTrajectory ?
                        updateParticle(particleStates[p], dtResamp, particleTrajectories.get(p))
                        : updateParticle(particleStates[p], dtResamp, null);

                double newWeight = Math.exp(newLogWeight);

                particleWeights[p] = newWeight;
                sumOfWeights += newWeight;
            }

            // Update marginal likelihood estimate
            thisLogP += Math.log(sumOfWeights / nParticles);

            if (!(sumOfWeights > 0.0)) {
                return Double.NEGATIVE_INFINITY;
            }

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
     * @param prevTreeEventIdx index of previous tree event
     * @param dtResamp length of interval
     * @param particleTrajectory if non-null, add particle states to this trajectory
     *
     * @return log conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
                                  int prevTreeEventIdx, double dtResamp,
                                  List<EpidemicState> particleTrajectory) {
        double conditionalLogP = 0;

        while (true) {

            model.calculatePropensities(particleState);

            double infectionProp = model.propensities[EpidemicEvent.INFECTION];
            double allowedInfectProp = infectionProp
                    *(1.0 - lineages * (lineages - 1) / particleState.I / (particleState.I + 1));
            double forbiddenInfectProp = infectionProp - allowedInfectProp;

            double allowedRecovProp, forbiddenRecovProp;
            if (particleState.I > lineages) {
                allowedRecovProp = model.propensities[EpidemicEvent.RECOVERY];
                forbiddenRecovProp = 0.0;
            } else {
                allowedRecovProp = 0.0;
                forbiddenRecovProp = model.propensities[EpidemicEvent.RECOVERY];
            }

            double allowedEventProp = allowedInfectProp + allowedRecovProp;

            // Determine length of proposed leap
            double tau = model.getTau(epsilon, particleState, allowedInfectProp, allowedRecovProp);


            double nextModelEventTime = model.getNextModelEventTime(particleState);
            double trueDt = Math.min(tau, Math.min(nextModelEventTime, finalTreeEvent.time) - particleState.time);
            conditionalLogP += -trueDt * (model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE]
                    + forbiddenInfectProp + forbiddenRecovProp);

            EpidemicEvent infectEvent = new EpidemicEvent();
            infectEvent.type = EpidemicEvent.INFECTION;
            infectEvent.multiplicity = (int)Randomizer.nextPoisson(trueDt*allowedInfectProp);

            EpidemicEvent recovEvent = new EpidemicEvent();
            recovEvent.type = EpidemicEvent.RECOVERY;
            recovEvent.multiplicity = (int)Math.round(Randomizer.nextPoisson(trueDt*allowedRecovProp));

            model.incrementState(particleState, infectEvent);
            model.incrementState(particleState, recovEvent);

            if (conditionalLogP == Double.NEGATIVE_INFINITY
                    || !particleState.isValid() || particleState.I < lineages)
                return Double.NEGATIVE_INFINITY;

            if (nextModelEventTime < finalTreeEvent.time && particleState.time + tau > nextModelEventTime) {
                particleState.time = nextModelEventTime;
                particleState.intervalIdx += 1;
            } else {
                if (particleState.time + tau > finalTreeEvent.time)
                    break;
                else
                    particleState.time += tau;
            }
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

}
