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
import beast.core.Operator;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Combines the TrajectoryOperator with one or more additional "
        + "operators (e.g. tree or treeOrigin operators).  This is necessary in"
        + "the exact sampler, as modifying the tree will inevitibly result in"
        + "a tree which is incompatible with the existing epidemic trajectory.")
public class CombinedSIRTrajectoryOperator extends SIRTrajectoryOperator {

    public Input<List<Operator>> operatorsInput = new Input<List<Operator>>(
            "operator",
            "Operator to combine with trajectory operator.",
            new ArrayList<Operator>());
    
    @Override
    public double proposal() {
        return super.proposal();
    }
    
}
