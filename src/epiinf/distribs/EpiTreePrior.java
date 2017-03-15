/*
 * Copyright (C) 2017 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.core.Input;
import beast.evolution.tree.TreeDistribution;
import epiinf.EpidemicTrajectory;
import epiinf.models.EpidemicModel;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EpiTreePrior extends TreeDistribution {
    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Input.Validate.REQUIRED);

    abstract public EpidemicTrajectory getConditionedTrajectory();
}
