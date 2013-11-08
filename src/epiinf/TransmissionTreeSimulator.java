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
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import com.google.common.collect.Lists;
import java.util.List;

/**
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
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
    
    public Input<Double> samplingProbInput = new Input<Double>(
            "samplingProbability",
            "Sampling probability: affects number of leaves in tree.",
            Validate.REQUIRED);

    private Tree tree;
    private RealParameter treeOrigin;
    private EpidemicTrajectory traj;
    private double samplingProb;
    
    public TransmissionTreeSimulator() { }
    
    @Override
    public void initAndValidate() throws Exception {
        tree = treeInput.get();
        treeOrigin = treeOriginInput.get();
        traj = trajInput.get();
        samplingProb = samplingProbInput.get();
    }
    
    @Override
    public void initStateNodes() throws Exception {
        Node transmissionRoot = new Node();
        
        List<Node> leaves = Lists.newArrayList();
        for (EpidemicEvent event : traj.getEventList()) {
            if (event.type == EpidemicEvent.EventType.RECOVERY
                    && Randomizer.nextDouble()<samplingProb) {
                Node leaf = new Node();
                leaf.setHeight(event.time);
                leaves.add(leaf);
            }
        }
        
        // Order leaves from most recent to oldest.
        leaves = Lists.reverse(leaves);
        
        
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(tree);
        stateNodes.add(treeOrigin);
    }
    
}
