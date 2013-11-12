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
import java.util.List;

/**
 * Simulate an epidemic trajectory under a stochastic SIR model.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate an epidemic trajectory under a stochastic SIR model.")
public class EpidemicTrajectorySimulator extends EpidemicTrajectory implements StateNodeInitialiser {
    
    public Input<Integer> S0Input = new Input<Integer>(
            "S0", "Initial susceptible count.",
            Validate.REQUIRED);

    public Input<Double> infectionRateInput = new Input<Double>(
            "infectionRate", "Rate of infection (per susceptible per infected)",
            Validate.REQUIRED);

    public Input<Double> recoveryRateInput = new Input<Double>(
            "recoveryRate", "Rate of recovery (per infected)",
            Validate.REQUIRED);
    
    double infectionRate, recoveryRate;
    
    public EpidemicTrajectorySimulator() { }
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        infectionRate = infectionRateInput.get();
        recoveryRate = recoveryRateInput.get();
        
        initialState = new EpidemicState(S0Input.get(), 1, 0);
        stateListDirty = true;
        
        simulate();
    }
    
    private void simulate() {
        eventList.clear();
        
        EpidemicState thisState = initialState.copy();
        double t = 0.0;

        while (thisState.I>0) {
            
            double infProp = thisState.S*thisState.I*infectionRate;
            double recProp = thisState.I*recoveryRate;
            double totalProp = infProp + recProp;

            t += Randomizer.nextExponential(totalProp);

            EpidemicEvent nextEvent = new EpidemicEvent();
            nextEvent.time = t;

            if (Randomizer.nextDouble()*totalProp < infProp) {
                nextEvent.type = EpidemicEvent.EventType.INFECTION;
                thisState.S -= 1;
                thisState.I += 1;
            } else {
                nextEvent.type = EpidemicEvent.EventType.RECOVERY;
                thisState.I -= 1;
                thisState.R += 1;
            }
            
            eventList.add(nextEvent);
        }
        
        stateListDirty = true;
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
