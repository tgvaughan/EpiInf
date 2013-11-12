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
import beast.core.parameter.RealParameter;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import epiinf.TreeEventList;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Exact probability of tree given ")
public class TreeDensity extends Distribution {
    
    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory object.", Validate.REQUIRED);
    
    EpidemicTrajectory trajectory;
    TreeEventList treeEventList;
    
    public TreeDensity() { }
    
    @Override
    public void initAndValidate() {
        trajectory = trajectoryInput.get();
        treeEventList = treeEventListInput.get();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
        List<TreeEventList.TreeEvent> revTreeEventList =
                Lists.reverse(treeEventList.getEventList());
        List<EpidemicEvent> revEventList = Lists.reverse(trajectory.getEventList());
        List<EpidemicState> revStateList = Lists.reverse(trajectory.getStateList());
        
        int idx = 0;

        int k = 0;
        
        for (TreeEventList.TreeEvent treeEvent : revTreeEventList) {
            
            while (idx<revEventList.size() &&
                    !TreeEventList.eventsMatch(treeEvent, revEventList.get(idx))) {
                
                // Incoporate probability of no effect on tree
                if (revEventList.get(idx).type == EpidemicEvent.EventType.INFECTION && k>1) {
                    int N = revStateList.get(idx).I;
                    logP += Math.log(1.0 - k*(k-1)/((double)N*(N-1)));
                }
                
                idx += 1;
            }
            
            // Check for incompatible tree and trajectory
            if (idx==revEventList.size()) {
                logP = Double.NEGATIVE_INFINITY;
                return logP;
            }
            
            // Incorporate probability of effect on tree
            if (treeEvent.type == TreeEventList.TreeEventType.SAMPLE)
                k += 1; // Sample
            else {
                int N = revStateList.get(idx).I;
                logP += Math.log(k*(k-1)/((double)N*(N-1)));
                k -= 1; // Coalescence
            }
            
        }
        
        return logP;
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

}
