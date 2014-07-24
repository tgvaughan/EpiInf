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
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.util.Randomizer;
import epiinf.models.EpidemicModel;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Simulate an epidemic trajectory under a stochastic SIR model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate an epidemic trajectory under a stochastic SIR model.")
public class TrajectorySimulator extends EpidemicTrajectory implements StateNodeInitialiser {
    
    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<Double> durationInput = new Input<>(
            "maxDuration", "Maximum duration of epidemic to simulate. "
                    + "Defaults to infinity.", Double.POSITIVE_INFINITY);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName",
            "Optional name of file to write simulated trajectory to.");
    
    public Input<Integer> nSerialSamplesInput = new Input<>(
            "nSerialSamples",
            "Optional number of serial samples.", 0);
    
    public Input<List<Double>> samplingTimesInput = new Input<>(
            "samplingTime", "Times (monotonically increasing) for contemporaneous sampling.",
            new ArrayList<>());
    
    public Input<List<Integer>> sampleSizesInput = new Input<>(
            "sampleSize", "Sizes of contemporaneous samples.",
            new ArrayList<>());
    
    EpidemicModel model;
    double duration;
    
    public TrajectorySimulator() { }
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        
        model = modelInput.get();
        duration = durationInput.get();
        
        if (samplingTimesInput.get().size() != sampleSizesInput.get().size())
            throw new IllegalArgumentException(
                    "Number of sample times and sample size entries do not match.");
        
        simulate();
        
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                    dumpTrajectory(ps);
            } catch (FileNotFoundException ex) {
                Logger.getLogger(TrajectorySimulator.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    
    private void simulate() {
        eventList.clear();
        stateList.clear();
        
        stateList.add(model.getInitialState());
        
        double t = 0.0;
        EpidemicState currentState = model.getInitialState();
        
        // Simulate portion of trajectory incorporating serial samples.
        for (int i=0; i<samplingTimesInput.get().size(); i++) {
            double sampleTime = samplingTimesInput.get().get(i);
            
            model.generateTrajectory(currentState, t, sampleTime);
            eventList.addAll(model.getEventList());
            stateList.addAll(model.getStateList());
            
            currentState = model.getStateList().get(model.getStateList().size()-1).copy();
            t = sampleTime;
            
            // Perform contemporaneous sampling
            for (int s=0; s<sampleSizesInput.get().get(i); s++) {
                model.incrementState(currentState, EpidemicEvent.Type.SAMPLE);
                EpidemicEvent samplingEvent = new EpidemicEvent();
                samplingEvent.time = t;
                samplingEvent.type = EpidemicEvent.Type.SAMPLE;
                eventList.add(samplingEvent);
                stateList.add(currentState.copy());
            }
        }
        
        // Simulate rest of trajectory
        model.generateTrajectory(currentState, t, duration);
        eventList.addAll(model.getEventList());
        stateList.addAll(model.getStateList());
        
        // Sample nSamples uniformly from amongst recovery events.
        if (nSerialSamplesInput.get()>0) {
            List<EpidemicEvent> recovEvents = new ArrayList<>();
            for (EpidemicEvent event : eventList)
                if (event.type == EpidemicEvent.Type.RECOVERY)
                    recovEvents.add(event);
            
            for (int i=0; i<nSerialSamplesInput.get(); i++) {
                EpidemicEvent event = recovEvents.get(Randomizer.nextInt(recovEvents.size()));
                event.type = EpidemicEvent.Type.SAMPLE;
                recovEvents.remove(event);
            }
        }
    }
    
    @Override
    public void initStateNodes() throws Exception {
        //simulate();
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
    
}
