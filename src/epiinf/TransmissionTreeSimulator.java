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

    public Input<EpidemicTrajectory> trajInput = new Input<>(
            "epidemicTrajectory",
            "Epidemic trajectory object.",
            Validate.REQUIRED);
    
    public Input<EpidemicModel> modelInput;
    
    public Input<String> fileNameInput = new Input<>(
            "fileName", "Name of file to save Newick representation of tree to.");
    
    public Input<Boolean> truncateTrajectoryInput = new Input<>(
            "truncateTrajectory",
            "Truncate trajectory at most recent sample. (Default true.)", true);

    private Tree tree;
    private EpidemicTrajectory traj;
    private boolean truncateTrajectory;
    
    public TransmissionTreeSimulator() {this.modelInput = new Input<>(
            "model",
            "Epidemic model.",
            Validate.REQUIRED);
 }
    
    @Override
    public void initAndValidate() throws Exception {
        tree = treeInput.get();
        traj = trajInput.get();
        
        truncateTrajectory = truncateTrajectoryInput.get();
    }
    
    @Override
    public void initStateNodes() throws Exception {

        int nSamples = 0;
        double youngestSamp = 0.0;
        for (EpidemicEvent event : traj.getEventList()) {
            if (event.type == EpidemicEvent.Type.SAMPLE) {
                nSamples += 1;
                youngestSamp = Math.max(event.time, youngestSamp);
            }
        }

        int nextLeafNr = 0;
        int nextInternalID = nSamples;
        
        List<Node> activeNodes = Lists.newArrayList();
        List<EpidemicEvent> revEventList = Lists.reverse(traj.getEventList());
        List<EpidemicState> revStateList = Lists.reverse(traj.getStateList());
        
        int samplesSeen = 0;
        for (int eidx=0; eidx<revEventList.size(); eidx++) {
            
            EpidemicEvent event = revEventList.get(eidx);
            EpidemicState state = revStateList.get(eidx);
            
            if (event.type == EpidemicEvent.Type.SAMPLE) {                
                Node leaf = new Node();
                leaf.setHeight(youngestSamp-event.time);
                leaf.setNr(nextLeafNr++);
                activeNodes.add(leaf);
                
                samplesSeen += 1;
            } else {
                
                if (event.type == EpidemicEvent.Type.INFECTION) {
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
            if (samplesSeen==nSamples && activeNodes.size()<2)
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
        
        // Write tree to disk if requested:
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {
                String newick = tree.toString().concat(";");
                ps.println(newick.replace(":0.0;", ":" + (youngestSamp-root.getHeight()) + ";"));
            }
        }
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(tree);
        stateNodes.add(traj);
    }
    
}
