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
import epiinf.EpidemicTrajectory;
import epiinf.TreeEventList;
import epiinf.models.EpidemicModel;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Uses SMC to propose a new trajectory conditional on the tree and "
        + "the epidemic origin.")
public class SMCTrajectoryOperator extends Operator {

    public Input<TreeEventList> treeEventListInput = new Input<TreeEventList>(
            "treeEventList", "Tree event list.", Input.Validate.REQUIRED);

    public Input<EpidemicTrajectory> trajInput = new Input<EpidemicTrajectory>(
            "epidemicTrajectory", "Epidemic trajectory.", Input.Validate.REQUIRED);
    
    public Input<EpidemicModel> modelInput = new Input<EpidemicModel>(
            "model", "Epidemic model.", Input.Validate.REQUIRED);
    
    @Override
    public double proposal() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
