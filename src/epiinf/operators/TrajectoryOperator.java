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
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.models.EpidemicModel;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which proposes new trajectories conditional on "
        + "the transmission tree and the tree origin.")
public class TrajectoryOperator extends Operator {
    
    public Input<TreeEventList> treeEventListInput = new Input<>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajInput = new Input<>(
            "epidemicTrajectory", "Epidemic trajectory.", Validate.REQUIRED);
    
    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    int count = 0;

    @Override
    public void initAndValidate() { }
    
    @Override
    public double proposal() {
        
        EpidemicModel model = modelInput.get();
        
        // Initialize HR with reverse move density
        double logHR = getReverseMoveDensity();
        
        EpidemicTrajectory traj = trajInput.get();
        
        List<EpidemicEvent> eventList = Lists.newArrayList();
        List<EpidemicState> stateList = Lists.newArrayList();
        List<TreeEvent> treeEventList = treeEventListInput.get().getEventList();

        EpidemicState thisState = model.getInitialState();
        stateList.add(thisState.copy());
        
        double lastTime = 0.0;
        for (TreeEvent treeEvent : treeEventList) {
            
            logHR -= model.generateTrajectory(thisState, lastTime, treeEvent.time);
            
            eventList.addAll(model.getEventList());
            stateList.addAll(model.getStateList());
            
            thisState = stateList.get(stateList.size()-1).copy();
            
            EpidemicEvent newEvent = new EpidemicEvent();
            newEvent.time = treeEvent.time;
            if (treeEvent.type == TreeEvent.Type.COALESCENCE) {
                newEvent.type = EpidemicEvent.Type.INFECTION;
                model.incrementState(thisState, newEvent.type);
            } else {
                newEvent.type = EpidemicEvent.Type.SAMPLE;
                model.incrementState(thisState, newEvent.type);
            }
            
            if (!thisState.isValid())
                return Double.NEGATIVE_INFINITY;
            
            eventList.add(newEvent);
            stateList.add(thisState.copy());
            
            lastTime = treeEvent.time;
        }
        
        traj.setEventAndStateList(eventList, stateList);
        
        return logHR;
    }
    
    /**
     * Obtain probability density of proposing existing trajectory given
     * current tree.
     * 
     * @return log probability density
     */
    private double getReverseMoveDensity() {
        
        EpidemicModel model = modelInput.get();
        
        List<TreeEvent> treeEventList = treeEventListInput.get().getEventList();
        List<EpidemicEvent> eventList = trajInput.get().getEventList();
        List<EpidemicState> stateList = trajInput.get().getStateList();
        
        double logP = 0.0;
        
        int lastIdx= -1;
        double lastTime = 0.0;
        for (TreeEvent treeEvent : treeEventList) {

            int nextIdx = lastIdx+1;
            while (nextIdx<eventList.size() &&
                    !model.eventsMatch(treeEvent, eventList.get(nextIdx)))
                nextIdx += 1;
            
            if (nextIdx>=eventList.size())
                return Double.NEGATIVE_INFINITY;
           
            int fromIdx;
            if (lastIdx>=0)
                fromIdx = lastIdx+1;
            else
                fromIdx = 0;
                
            logP += model.getPathProbability(lastTime, treeEvent.time,
                        stateList.get(fromIdx),
                        eventList.subList(fromIdx, nextIdx));
            
            lastTime = treeEvent.time;
            lastIdx = nextIdx;
        }
        
        return logP;
    }
    
}
