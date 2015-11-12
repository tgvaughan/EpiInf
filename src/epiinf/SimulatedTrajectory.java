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

package epiinf;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.util.Randomizer;
import epiinf.models.EpidemicModel;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simulate an epidemic trajectory under a stochastic SIR model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate an epidemic trajectory under a stochastic SIR model.")
public class SimulatedTrajectory extends EpidemicTrajectory {
    
    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<Double> durationInput = new Input<>(
            "maxDuration", "Maximum duration of epidemic to simulate. "
                    + "Defaults to infinity.", Double.POSITIVE_INFINITY);

    public Input<Function> originInput = new Input<>(
            "origin", "Origin with respect to most recent sample in tree. " +
            "If provided, trimes will be logged as ages before most recent " +
            "sample.");

    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Optional name of file to write simulated trajectory to.");
    
    EpidemicModel model;
    double duration;
    
    public SimulatedTrajectory() { }
    
    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();
        
        model = modelInput.get();
        duration = durationInput.get();

        simulate();
        
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                    dumpTrajectory(ps);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(SimulatedTrajectory.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }


    /**
     * Simulate an epidemic using the provided model.
     *
     */
    public void simulate() {

        double endTime;
        if (origin != null)
            endTime = origin;
        else
            endTime = duration;

        eventList.clear();
        stateList.clear();

        EpidemicState thisState = model.getInitialState();
        stateList.add(model.getInitialState());

        thisState.time = 0;

        while (true) {
            model.calculatePropensities(thisState);

            double totalPropensity = model.propensities[EpidemicEvent.INFECTION]
                    + model.propensities[EpidemicEvent.RECOVERY]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_REMOVE]
                    + model.propensities[EpidemicEvent.PSI_SAMPLE_NOREMOVE];

            double dt;
            if (totalPropensity>0.0)
                dt = Randomizer.nextExponential(totalPropensity);
            else
                dt = Double.POSITIVE_INFINITY;

            thisState.time += dt;

            EpidemicEvent nextEvent = new EpidemicEvent();

            double nextModelEventTime = model.getNextModelEventTime(thisState);
            if (nextModelEventTime < endTime && thisState.time > nextModelEventTime) {
                ModelEvent event = model.getNextModelEvent(thisState);

                if (event.type == ModelEvent.Type.RHO_SAMPLING) {
                    nextEvent.type = EpidemicEvent.RHO_SAMPLE;

                    // Got to be a better way of sampling from a binomial distribution
                    nextEvent.multiplicity = 0;
                    for (int i = 0; i < thisState.I; i++) {
                        if (Randomizer.nextDouble() < event.rho)
                            nextEvent.multiplicity += 1;
                    }


                    nextEvent.time = event.time;
                    thisState.time = event.time;

                    model.incrementState(thisState, nextEvent);
                    eventList.add(nextEvent);
                    stateList.add(thisState.copy());
                }

                thisState.intervalIdx += 1;

                continue;
            }

            if (thisState.time>=endTime)
                break;


            nextEvent.time = thisState.time;

            double u = totalPropensity*Randomizer.nextDouble();

            for (int type = 0; type<EpidemicEvent.nTypes; type++) {
                u -= model.propensities[type];

                if (u<0) {
                    nextEvent.type = type;
                    break;
                }
            }

            model.incrementState(thisState, nextEvent);
            eventList.add(nextEvent);
            stateList.add(thisState.copy());
        }
    }
    
    @Override
    public void log(int nSample, PrintStream out) {
        if (originInput.get() != null)
            origin = originInput.get().getArrayValue();

        simulate();

        super.log(nSample, out);
    }
}
