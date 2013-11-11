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

package epiinf.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.TreeEventList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which proposes new trajectories conditional on "
        + "the transmission tree and the tree origin.")
public class TrajectoryOperator extends Operator {
    
    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate", "Rate of infection.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "recoveryRate", "Rate of recovery.", Validate.REQUIRED);
    
    public Input<IntegerParameter> S0Input = new Input<IntegerParameter>(
            "S0", "Initial number of susceptibles.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory.", Validate.REQUIRED);

    
    @Override
    public double proposal() {
        double logHR = 0.0;
        
        EpidemicTrajectory traj = trajInput.get();
        double infectionRate = infectionRateInput.get().getValue();
        double recoveryRate = recoveryRateInput.get().getValue();
        
        List<EpidemicEvent> eventList = Lists.newArrayList();
        List<TreeEventList.TreeEvent> treeEventList = treeEventListInput.get().getEventList();

        EpidemicState thisState = new EpidemicState(S0Input.get().getValue(), 1, 0);
        traj.setInitialState(thisState);
        
        double t = 0.0;
        for (TreeEventList.TreeEvent treeEvent : treeEventList) {
            while (true) {
                double infProp = infectionRate*thisState.S*thisState.I;
                double recProp = recoveryRate*thisState.I;
                double totalProp = infProp + recProp;
                
                t += Randomizer.nextExponential(totalProp);
                if (t > treeEvent.time) {
                    t = treeEvent.time;
                    break;
                }
                
                EpidemicEvent newEvent = new EpidemicEvent();
                newEvent.time = t;
                
                if (Randomizer.nextDouble()*totalProp < infProp) {
                    // Infection
                    newEvent.type = EpidemicEvent.EventType.INFECTION;
                    eventList.add(newEvent);
                    
                    thisState.S -= 1;
                    thisState.I += 1;
                } else {
                    // Recovery
                    newEvent.type = EpidemicEvent.EventType.RECOVERY;
                    eventList.add(newEvent);
                    
                    thisState.I -= 1;
                    thisState.R += 1;
                }
            }
            
            EpidemicEvent newEvent = new EpidemicEvent();
            newEvent.time = treeEvent.time;
            if (treeEvent.type == TreeEventList.TreeEventType.COALESCENCE)
                newEvent.type = EpidemicEvent.EventType.INFECTION;
            else
                newEvent.type = EpidemicEvent.EventType.RECOVERY;
        }
        
        return logHR;
    }
    
}
