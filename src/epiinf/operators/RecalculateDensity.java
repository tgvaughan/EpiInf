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

package epiinf.operators;

import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;

import java.io.PrintStream;

/**
 * Dummy Operator that causes the tree prior to be recomputed.  Probability of selecting this
 * operator must go to zero for the sampler to remain exact.
 *
 * Created by Tim Vaughan <tgvaughan@gmail.com> on 24/04/17.
 */
public class RecalculateDensity extends Operator implements Loggable {

    public Input<StateNode> stateNodeInput = new Input("stateNode", "StateNode to reevaluate density of.", Input.Validate.REQUIRED);
    public Input<Double> iterationsInput = new Input<>("iterations", "Number of iterations this operator should be used for.");

    boolean stillOn;

    public RecalculateDensity() { }

    @Override
    public void initAndValidate() {
        stillOn = true;
    }

    @Override
    public double proposal() {
        if (!stillOn)
            return Double.NEGATIVE_INFINITY;

        stateNodeInput.get().startEditing(this);
        return Double.POSITIVE_INFINITY;
    }

    @Override
    public void init(PrintStream out) {
    }

    @Override
    public void log(long sample, PrintStream out) {
        if (iterationsInput.get() != null)
            stillOn = sample < iterationsInput.get();
    }

    @Override
    public void close(PrintStream out) {
    }
}
