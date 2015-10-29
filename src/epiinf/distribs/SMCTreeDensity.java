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

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.math.Binomial;
import beast.math.GammaFunction;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.*;
import epiinf.models.EpidemicModel;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
    + "parameters.")
public class SMCTreeDensity extends Distribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);

    public Input<TreeEventList> treeEventListInput = new Input<>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);

    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    EpidemicModel model;
    TreeEventList treeEventList;
    int nParticles;

    // DEBUG
//    PrintStream debugOut;

    public SMCTreeDensity() {
    }

    @Override
    public void initAndValidate() throws Exception {
        model = modelInput.get();
        treeEventList = treeEventListInput.get();
        nParticles = nParticlesInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {

        logP = 0.0;

        if (treeEventList.getEventList().get(0).time < 0) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }

        List<ModelEvent> modelEventList = model.getModelEventList();

        List<Double> particleWeights = Lists.newArrayList();
        List<EpidemicState> particleStates = Lists.newArrayList();
        List<EpidemicState> particleStatesNew = Lists.newArrayList();

        // Initialize particles
        for (int p = 0; p < nParticles; p++)
            particleStates.add(model.getInitialState());

        double t = 0.0;
        int k = 1;
        for (TreeEvent treeEvent : treeEventList.getEventList()) {

            // Update particles
            particleWeights.clear();
            double sumOfWeights = 0.0;
            for (int p = 0; p < nParticles; p++) {

                double newWeight = updateParticle(particleStates.get(p), t, k, treeEvent, modelEventList);

                particleWeights.add(newWeight);
                sumOfWeights += newWeight;
            }

            // Update marginal likelihood estimate
            logP += Math.log(sumOfWeights / nParticles);

            if (!(sumOfWeights > 0.0))
                return Double.NEGATIVE_INFINITY;

            // Sample particle with replacement
            particleStatesNew.clear();
            for (int p = 0; p < nParticles; p++) {
                double u = Randomizer.nextDouble() * sumOfWeights;

                int pChoice;
                for (pChoice = 0; pChoice < nParticles; pChoice++) {
                    u -= particleWeights.get(pChoice);
                    if (u < 0.0)
                        break;
                }

                if (pChoice == nParticles)
                    System.err.println("sumOfWeights: " + sumOfWeights);

                particleStatesNew.add(particleStates.get(pChoice).copy());
            }

            // Switch particleStates and particleStatesNew
            List<EpidemicState> temp = particleStates;
            particleStates = particleStatesNew;
            particleStatesNew = temp;

            // Update lineage counter
            if (treeEvent.type == TreeEvent.Type.COALESCENCE)
                k += 1;
            else
                k -= treeEvent.multiplicity;

            // Update start interval time
            t = treeEvent.time;
        }

        return logP;
    }

    /**
     * Updates weight and state of particle.
     *
     * @param particleState
     * @param startTime
     * @param lineages
     * @param finalTreeEvent
     * @return conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
                                  double startTime, int lineages, TreeEvent finalTreeEvent,
                                  List<ModelEvent> modelEventList) {
        double conditionalP = 1.0;

        double t = startTime;

        // Trim off any model events which occurred before the start of the tree
        while (!modelEventList.isEmpty() && modelEventList.get(0).time<startTime) {
            if (modelEventList.get(0).type == ModelEvent.Type.RHO_SAMPLING)
                return 0.0;

            modelEventList.remove(0);
        }

        double nextModelEventTime = Double.POSITIVE_INFINITY;
        if (!modelEventList.isEmpty())
            nextModelEventTime = modelEventList.get(0).time;

        while (true) {
            model.calculatePropensities(particleState);

            double infectionProp = model.propensities.get(EpidemicEvent.Type.INFECTION);
            double allowedRecovProp, forbiddenRecovProp;
            if (particleState.I > lineages) {
                allowedRecovProp = model.propensities.get(EpidemicEvent.Type.RECOVERY);
                forbiddenRecovProp = 0.0;
            } else {
                allowedRecovProp = 0.0;
                forbiddenRecovProp = model.propensities.get(EpidemicEvent.Type.RECOVERY);
            }
            double allowedEventProp = infectionProp + allowedRecovProp;

            // Determine size of time increment

            double dt;
            if (allowedEventProp > 0.0)
                dt = Randomizer.nextExponential(allowedEventProp);
            else
                dt = Double.POSITIVE_INFINITY;

            // Condition against psi-sampling and illegal recovery within interval
            double trueDt = Math.min(dt, Math.min(nextModelEventTime, finalTreeEvent.time) - t);
            conditionalP *= Math.exp(-trueDt*(model.propensities.get(EpidemicEvent.Type.PSI_SAMPLE) + forbiddenRecovProp));

            // Increment time
            t += dt;

            // Deal with model events (rho sampling and rate shifts)
            if (nextModelEventTime < finalTreeEvent.time && t > nextModelEventTime) {
                if (modelEventList.get(0).type == ModelEvent.Type.RHO_SAMPLING)
                    return 0.0;
                else {
                    t = nextModelEventTime;
                    modelEventList.remove(0);
                    if (modelEventList.isEmpty())
                        nextModelEventTime = Double.POSITIVE_INFINITY;
                    else
                        nextModelEventTime = modelEventList.get(0).time;

                    continue;
                }
            }

            // Stop here if we're past the end of the tree interval
            if (t > finalTreeEvent.time)
                break;

            EpidemicEvent event = new EpidemicEvent();
            if (allowedEventProp*Randomizer.nextDouble() < infectionProp) {
                event.type = EpidemicEvent.Type.INFECTION;
            } else
                event.type = EpidemicEvent.Type.RECOVERY;


            // Condition against infection events that produce coalescences not
            // observed in tree.
            if (event.type == EpidemicEvent.Type.INFECTION)
                conditionalP *= 1.0 - lineages * (lineages - 1) / particleState.I / (particleState.I + 1);

            model.incrementState(particleState, event);

            // Early exit if invalid state:
            if (conditionalP == 0)
                return 0.0;

        }

        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.calculatePropensities(particleState);
            model.incrementState(particleState, EpidemicEvent.Infection);
            conditionalP *= 2.0 / particleState.I / (particleState.I - 1)
                    * model.propensities.get(EpidemicEvent.Type.INFECTION);
        } else {

            double sampleProb;
            if (model.timesEqual(finalTreeEvent.time, nextModelEventTime) && modelEventList.get(0).type == ModelEvent.Type.RHO_SAMPLING) {

                sampleProb = 0.0;

                // If model contains a rho sampling event at this time, calculate the probability
                // of sampling the number of samples in finalTreeEvent given the current
                // state.
                for (int i = 0; i < model.rhoSamplingProbInput.get().getDimension(); i++) {
                    double rhoProb = model.rhoSamplingProbInput.get().getValue(i);
                    double rhoTime = model.rhoSamplingTimeInput.get().getValue(i);

                    if (Math.abs(rhoTime - finalTreeEvent.time) < model.getTolerance()) {
                        int I = (int) Math.round(particleState.I);
                        int k = finalTreeEvent.multiplicity;
                        sampleProb += Binomial.choose(I, k)
                                * Math.pow(rhoProb, k) * Math.pow(1.0 - rhoProb, I - k);
                    }
                }

                model.incrementState(particleState,
                        EpidemicEvent.MultipleRhoSamples(finalTreeEvent.multiplicity));

            } else {
                if (model.psiSamplingRateInput.get() != null && finalTreeEvent.multiplicity == 1) {
                    model.calculatePropensities(particleState);
                    sampleProb = model.propensities.get(EpidemicEvent.Type.PSI_SAMPLE);
                    model.incrementState(particleState, EpidemicEvent.PsiSample);
                } else {
                    // No explicit sampling process
                    sampleProb = 1.0;
                    model.incrementState(particleState,
                            EpidemicEvent.MultipleOtherSamples(finalTreeEvent.multiplicity));
                }

            }

            conditionalP *= sampleProb * Math.exp(GammaFunction.lnGamma(1 + finalTreeEvent.multiplicity));
        }

        if (!particleState.isValid())
            return 0.0;
        else
            return conditionalP;
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
    public boolean isStochastic() {
        return true;
    }
}
