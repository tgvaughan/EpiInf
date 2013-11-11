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
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which proposes new trajectories conditional on "
        + "the transmission tree and the tree origin.")
public class TrajectoryOperator extends Operator {
    
    public Input<Tree> treeInput = new Input<Tree>(
            "tree", "Transmission tree.", Validate.REQUIRED);
    
    public Input<RealParameter> treeOriginInput = new Input<RealParameter>(
            "treeOrigin", "Difference between time of MRCA and start of "
                    + "epidemic.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate", "Rate of infection.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "recoveryRate", "Rate of recovery.", Validate.REQUIRED);

    @Override
    public double proposal() {
        double logHR = 0.0;
        
        return logHR;
    }
    
}
