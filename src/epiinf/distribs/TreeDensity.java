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
import beast.evolution.tree.Tree;
import epiinf.EpidemicTrajectory;
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

    Tree tree;
    EpidemicTrajectory trajectory;
    
    public TreeDensity() { }
    
    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        trajectory = trajectoryInput.get();
    }
    
    @Override
    public double calculateLogP() {
        logP = 0.0;
        
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
