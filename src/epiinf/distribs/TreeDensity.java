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
import epiinf.EpidemicTrajectory;
import java.util.ArrayList;
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
    
    enum TreeEventType { SAMPLE, COALESCENCE };
    class TreeEvent {
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
        
        
        return logP;
    }
    
    /**
     * Ensure list of tree events is up to date.
     */
    private void updateTreeEventList() {
        treeEventList.clear();

        for (Node node : tree.getNodesAsArray()) {
            TreeEvent event = new TreeEvent();
            if (node.isLeaf())
                event.type = TreeEventType.SAMPLE;
            else
                event.type = TreeEventType.COALESCENCE;
            
        }
        
    }
    
    private double getTimeFromHeight(double height) {
        
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
