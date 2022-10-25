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

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Operator which acts on a RealParameter, maintaining the " +
        "order of all elements")
public class ChangeTimesOperator extends Operator {

    public Input<RealParameter> changeTimesInput = new Input<>(
            "changeTimes", "Parameter describing change times.",
            Input.Validate.REQUIRED);

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor", "Determines range of scaling operation.",
            0.8);

    @Override
    public double proposal() {

        RealParameter changeTimes = changeTimesInput.get();

        double minf = Math.min(scaleFactorInput.get(), 1.0/scaleFactorInput.get());
        double f = minf + (1.0/minf - minf)* Randomizer.nextDouble();

        // Select element:

        int i = Randomizer.nextInt(changeTimes.getDimension());
        double newValue = f*changeTimes.getValue(i);

        if (newValue < changeTimes.getLower()
                || newValue > changeTimes.getUpper()
                || i>0 && newValue < changeTimes.getValue(i-1)
                || i<changeTimes.getDimension()-1 && newValue > changeTimes.getValue(i+1))
            return Double.NEGATIVE_INFINITY;

        changeTimes.setValue(i, newValue);

        return -Math.log(f);
    }

    @Override
    public void initAndValidate() { }
}
