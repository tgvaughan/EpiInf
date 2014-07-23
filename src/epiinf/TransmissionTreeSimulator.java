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

import epiinf.models.EpidemicModel;
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
import java.io.PrintStream;
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

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree to initialize",
            Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<>(
            "treeOrigin",
            "Difference between start of epidemic and MRCA.",
            Validate.REQUIRED);
    
    public Input<EpidemicTrajectory> trajInput = new Input<>(
            "epidemicTrajectory",
            "Epidemic trajectory object.",
            Validate.REQUIRED);
    
    public Input<EpidemicModel> modelInput;
    
    public Input<Integer> nLeavesInput = new Input<>(
            "nLeaves",
            "Fixed number of leaves in tree.");
    
    public Input<Double> sampleBeforeTimeInput = new Input<>(
            "sampleBeforeTime",
            "Sample leaves before this time, specified as a fraction of the "
                    + "total epidemic time.  Default 1.0.", 1.0);
    
    public Input<String> fileNameInput = new Input<>(
            "fileName", "Name of file to save Newick representation of tree to.");
    
    public Input<Boolean> truncateTrajectoryInput = new Input<>(
            "truncateTrajectory",
            "Truncate trajectory at most recent sample. (Default true.)", true);
    
    public Input<Boolean> contempSamplingInput = new Input<>(
            "contempSampling",
            "Leaves are sampled at the same time (end of trajectory or"
                    + " fraction given by sampleBeforeTime).", false);

    private Tree tree;
    private RealParameter treeOrigin;
    private EpidemicTrajectory traj;
    private EpidemicModel model;
    private int nLeaves = -1;
    private double sampleBeforeTime;
    private boolean truncateTrajectory;
    private boolean contempSampling;
    
    public TransmissionTreeSimulator() {this.modelInput = new Input<>(
            "model",
            "Epidemic model.",
            Validate.REQUIRED);
 }
    
    @Override
    public void initAndValidate() throws Exception {
        tree = treeInput.get();
        treeOrigin = treeOriginInput.get();
        traj = trajInput.get();
        model = modelInput.get();
        
        if (nLeavesInput.get() != null)
            nLeaves = nLeavesInput.get();
        
        sampleBeforeTime = sampleBeforeTimeInput.get();
        truncateTrajectory = truncateTrajectoryInput.get();
        contempSampling = contempSamplingInput.get();
    }
    
    @Override
    public void initStateNodes() throws Exception {
        
        List<EpidemicEvent> leafEvents = Lists.newArrayList();
        double epidemicDuration = traj.getEventList()
                .get(traj.getEventList().size()-1).time;
        
        double samplingTime = sampleBeforeTime*epidemicDuration;

        // Extract existing sampling events from tree
        List<EpidemicEvent> samplingEvents = Lists.newArrayList();
        for (EpidemicEvent event : traj.getEventList()) {
            if (event.type == model.getLeafEventType()
                    && event.time<=samplingTime)
                samplingEvents.add(event);
        }
        
        if (nLeaves>0) {
            // Check for invalid leaf count request.
            if (samplingEvents.size()<nLeaves)
                throw new IllegalArgumentException("Trajectory has fewer recovery"
                        + " events than value of nLeaves.");
            
            for (int i=0; i<nLeaves; i++) {
                int idx = Randomizer.nextInt(samplingEvents.size());
                leafEvents.add(samplingEvents.get(idx));
                samplingEvents.remove(idx);
            }
            
        } else {

            for (EpidemicEvent event : samplingEvents) {
                if (Randomizer.nextDouble()<model.getProbLeaf())
                    leafEvents.add(event);
            }
            
        }
        
        int nextLeafNr = 0;
        int nextInternalID = leafEvents.size();
        
        // Order leaf events from youngest to oldest.
        Collections.sort(leafEvents, (EpidemicEvent o1, EpidemicEvent o2) -> {
            if (o1.time > o2.time)
                return -1;
            
            if (o1.time < o2.time)
                return 1;
            
            return 0;
        });
        
        // Record time of youngest sample:
        double youngestSamp = leafEvents.get(0).time;
        
        List<Node> activeNodes = Lists.newArrayList();
        List<EpidemicEvent> revEventList = Lists.reverse(traj.getEventList());
        List<EpidemicState> revStateList = Lists.reverse(traj.getStateList());
        
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
                
                if (event.type == model.getCoalescenceEventType()) {
                    int k = activeNodes.size();
                    double N = state.I;
                    
                    double pCoalesce = k*(k-1)/(N*(N-1));
                    
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
        
        // Truncate trajectory at most recent sample if requested:
        if (truncateTrajectory) {
            while (revEventList.get(0).time>youngestSamp) {
                revEventList.remove(0);
                revStateList.remove(0);
            }
        }

        // Initialise state nodes
        Node root = activeNodes.get(0);
        tree.assignFromWithoutID(new Tree(root));
        treeOrigin.assignFromWithoutID(new RealParameter(String.valueOf(youngestSamp-root.getHeight())));
        
        // Write tree to disk if requested:
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                String newick = tree.toString().concat(";");
                ps.println(newick.replace(":0.0;", ":" + treeOrigin.getValue() + ";"));
            }
        }
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(tree);
        stateNodes.add(treeOrigin);
        stateNodes.add(traj);
    }
    
}
