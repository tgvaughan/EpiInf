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

import beast.core.Function;
import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.TreeDistribution;
import epiinf.EpidemicTrajectory;
import epiinf.ObservedEventsList;
import epiinf.models.EpidemicModel;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public abstract class EpiTreePrior extends TreeDistribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Input.Validate.REQUIRED);

    public Input<Function> finalSampleOffsetInput = new Input<>(
            "finalSampleOffset",
            "Difference in time between final sample and end of observation process", Input.Validate.REQUIRED);

    protected EpidemicModel model;
    protected ObservedEventsList observedEventsList;

    abstract public EpidemicTrajectory getConditionedTrajectory();

    /**
     * @return Epidemic model
     */
    public EpidemicModel getModel() {
        return modelInput.get();
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

    @Override
    protected boolean requiresRecalculation() {
        observedEventsList.makeDirty();
        return true;
    }

    @Override
    public void restore() {
        observedEventsList.makeDirty();
        super.restore();
    }

    @Override
    public boolean isStochastic() {
        return true;
    }

}
