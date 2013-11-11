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

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * State node initializer which simulates a transmission tree compatible with
 * the given epidemic trajectory.  The tree origin parameter is also initialized.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Simulate transmission tree conditional on epidemic trajectory.")
public class TransmissionTreeSimulator extends BEASTObject implements StateNodeInitialiser {

    public Input<Tree> treeInput = new Input<Tree>(
            "tree",
            "Tree to initialize",
            Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<RealParameter>(
            "treeOrigin",
            "Difference between start of epidemic and MRCA.",
            Validate.REQUIRED);
    
    public Input<EpidemicTrajectory> trajInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory",
            "Epidemic trajectory object.",
            Validate.REQUIRED);
    
    public Input<Integer> nLeavesInput = new Input<Integer>(
            "nLeaves",
            "Number of leaves in tree.",
            Validate.REQUIRED);

    private Tree tree;
    private RealParameter treeOrigin;
    private EpidemicTrajectory traj;
    private int nLeaves;
    
    public TransmissionTreeSimulator() { }
    
    @Override
    public void initAndValidate() throws Exception {
        tree = treeInput.get();
        treeOrigin = treeOriginInput.get();
        traj = trajInput.get();
        nLeaves = nLeavesInput.get();
    }
    
    @Override
    public void initStateNodes() throws Exception {
        
        List<EpidemicEvent> leafEvents = Lists.newArrayList();
        List<EpidemicEvent> recoveryEvents = Lists.newArrayList();
        
        for (EpidemicEvent event : traj.getEventList()) {
            if (event.type == EpidemicEvent.EventType.RECOVERY)
                recoveryEvents.add(event);
        }

        // Check for invalid leaf count request.
        if (recoveryEvents.size()<nLeaves)
            throw new IllegalArgumentException("Trajectory has fewer recovery"
                    + "events than value of nLeaves.");
        
        for (int i=0; i<nLeaves; i++) {
            int idx = Randomizer.nextInt(recoveryEvents.size());
            leafEvents.add(recoveryEvents.get(idx));
            recoveryEvents.remove(idx);
        }
        
        int nextLeafNr = 0;
        int nextInternalID = leafEvents.size();
        
        // Order leaf events from youngest to oldest.
        Collections.sort(leafEvents, new Comparator<EpidemicEvent>() {

            @Override
            public int compare(EpidemicEvent o1, EpidemicEvent o2) {
                if (o1.time > o2.time)
                    return -1;
                
                if (o1.time < o2.time)
                    return 1;
                
                return 0;
            }
        });
        
        // Record time of youngest sample:
        double youngestSamp = leafEvents.get(0).time;
        
        List<Node> activeNodes = Lists.newArrayList();
        List<EpidemicEvent> revEventList = Lists.reverse(traj.getEventList());
        List<EpidemicState> revStateList = Lists.reverse(traj.getStateList());
        
        //for (EpidemicEvent event : Lists.reverse(traj.getEventList())) {
        for (int eidx=0; eidx<revEventList.size(); eidx++) {
            
            EpidemicEvent event = revEventList.get(eidx);
            EpidemicState state = revStateList.get(eidx);
            
            if (!leafEvents.isEmpty() && leafEvents.get(0).equals(event)) {
                
                Node leaf = new Node();
                leaf.setHeight(youngestSamp-event.time);
                leaf.setNr(nextLeafNr++);
                activeNodes.add(leaf);
                leafEvents.remove(0);
                
            } else {
                
                if (event.type == EpidemicEvent.EventType.INFECTION) {
                    int k = activeNodes.size();
                    int N = state.I;
                    
                    double pCoalesce = k*(k-1)/((double)N*(N-1));
                    
                    if (Randomizer.nextDouble()<pCoalesce) {
                        int childIdx = Randomizer.nextInt(k);
                        Node child1 = activeNodes.get(childIdx);
                        activeNodes.remove(childIdx);
                        
                        childIdx = Randomizer.nextInt(k-1);
                        Node child2 = activeNodes.get(childIdx);
                        activeNodes.remove(childIdx);
                        
                        Node parent = new Node();
                        parent.addChild(child1);
                        parent.addChild(child2);
                        parent.setHeight(youngestSamp-event.time);
                        parent.setNr(nextInternalID++);
                        activeNodes.add(parent);
                    }
                }
            }
            
            // Stop when we reach the MRCA of all sampled events.
            if (leafEvents.isEmpty() && activeNodes.size()<2)
                break;
        }

        // Initialise state nodes
        Node root = activeNodes.get(0);
        tree.assignFromWithoutID(new Tree(root));
        treeOrigin.assignFromWithoutID(new RealParameter(String.valueOf(youngestSamp-root.getHeight())));
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(tree);
        stateNodes.add(treeOrigin);
    }
    
}
