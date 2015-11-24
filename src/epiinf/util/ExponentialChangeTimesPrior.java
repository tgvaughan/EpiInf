/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

package epiinf.util;

import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import org.apache.commons.math.distribution.Distribution;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ExponentialChangeTimesPrior extends beast.core.Distribution {

    public Input<RealParameter> changeTimesInput = new Input<>("changeTimes",
            "Change times parameter.", Input.Validate.REQUIRED);

    public Input<RealParameter> meanIntervalInput = new Input<>("meanIntervalLength",
            "Intervals lengths are exponentially distributed with this mean.",
            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
    }

    @Override
    public double calculateLogP() throws Exception {
        RealParameter changeTimes = changeTimesInput.get();
        double meanInterval = meanIntervalInput.get().getValue();

        logP = 0.0;
        double lastVal = 0.0;
        for (int i=0; i<changeTimes.getDimension(); i++) {
            double thisVal = changeTimes.getValue(i);
            logP += -(thisVal-lastVal)/meanInterval;
            lastVal = thisVal;
        }

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
        throw new UnsupportedOperationException("Sampling not supported.");
    }
}
