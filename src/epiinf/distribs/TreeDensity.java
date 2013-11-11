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
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.EpidemicTrajectory;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Exact probability of tree given ")
public class TreeDensity extends Distribution {
    
    public Input<Tree> treeInput = new Input<Tree>(
            "tree", "Tree to calculate the density of.", Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajectoryInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory object.", Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<RealParameter>(
            "treeOrigin", "Time between the epidemic start and the tree MRCA.",
            Validate.REQUIRED);
    
    Tree tree;
    EpidemicTrajectory trajectory;
    RealParameter treeOrigin;
    
    /**
     * Tolerance for deviation between tree node ages and trajectory events.
     */
    public static final double tolerance = 1e-10;
   
    public enum TreeEventType { SAMPLE, COALESCENCE }
    public class TreeEvent {
        Node node;
        TreeEventType type;
        double time;    
    }
    
    List<TreeEvent> treeEventList;
    
    public TreeDensity() { }
    
    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        trajectory = trajectoryInput.get();
        treeOrigin = treeOriginInput.get();
        
        treeEventList = new ArrayList<TreeEvent>();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
        updateTreeEventList();
        
        //List<TreeEvent> revTreeEventList = Lists.reverse(treeEventList);
        List<EpidemicEvent> revEventList = Lists.reverse(trajectory.getEventList());
        List<EpidemicState> revStateList = Lists.reverse(trajectory.getStateList());
        
        int idx = 0;

        int k = 0;
        
        for (TreeEvent treeEvent : treeEventList) {
            
            while (idx<revEventList.size()
                    && !eventsMatch(treeEvent, revEventList.get(idx))) {
                
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
            if (treeEvent.type == TreeEventType.SAMPLE)
                k += 1; // Sample
            else {
                int N = revStateList.get(idx).I;
                logP += Math.log(k*(k-1)/((double)N*(N-1)));
                k -= 1; // Coalescence
            }
            
        }
        
        return logP;
    }

    /**
     * Determine whether given tree event and epidemic event are compatible.
     * 
     * @param treeEvent
     * @param epiEvent
     * 
     * @return true if events are compatible
     */
    private boolean eventsMatch(TreeEvent treeEvent, EpidemicEvent epiEvent) {
        
        if (Math.abs(treeEvent.time-epiEvent.time)>TreeDensity.tolerance)
            return false;
        
        if ((treeEvent.type == TreeEventType.COALESCENCE)
                && (epiEvent.type != EpidemicEvent.EventType.INFECTION))
            return false;
        
        if ((treeEvent.type == TreeEventType.SAMPLE)
                && (epiEvent.type != EpidemicEvent.EventType.RECOVERY))
            return false;
        
        return true;
    }
    
    /**
     * Ensure list of tree events is up to date.
     */
    private void updateTreeEventList() {
        treeEventList.clear();

        // Assemble event list
        for (Node node : tree.getNodesAsArray()) {
            TreeEvent event = new TreeEvent();
            if (node.isLeaf())
                event.type = TreeEventType.SAMPLE;
            else
                event.type = TreeEventType.COALESCENCE;
       
            event.time = getTimeFromHeight(node.getHeight());
            event.node = node;
            
            treeEventList.add(event);
        }
        
        // Sort events in time-reversed order
        Collections.sort(treeEventList, new Comparator<TreeEvent>() {

            @Override
            public int compare(TreeEvent e1, TreeEvent e2) {
                if (e1.time < e2.time)
                    return 1;
                
                if (e1.time > e2.time)
                    return -1;
                
                return 0;
            }
            
        });
        
    }
    
    /**
     * Obtain absolute epidemic time corresponding to height on tree.
     * 
     * @param height
     * @return time
     */
    private double getTimeFromHeight (double height) {
        return treeOrigin.getValue() + tree.getRoot().getHeight() - height;
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
