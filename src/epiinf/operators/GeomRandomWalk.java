/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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
import beast.core.parameter.IntegerParameter;
import beast.util.Randomizer;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Proposes new value for integer parameter from geometric" +
        " distribution centred on current value.")
public class GeomRandomWalk extends Operator {

    public Input<IntegerParameter> parameterInput = new Input<>("parameter",
            "Integer parameter to scale.", Input.Validate.REQUIRED);

    IntegerParameter parameter;

    @Override
    public void initAndValidate() {
        parameter = parameterInput.get();
    }

    @Override
    public double proposal() {

        int dim = parameter.getDimension();
        int idx = dim > 1 ? Randomizer.nextInt(dim) : 0;

        int n = parameter.getValue(idx);

        if (n <=0)
            return Double.NEGATIVE_INFINITY;

        double p = 1/(double)n;

        int nPrime = (int)Randomizer.nextGeometric(p);

        if (nPrime <= 0
                || nPrime < parameterInput.get().getLower()
                || nPrime > parameterInput.get().getUpper())
            return Double.NEGATIVE_INFINITY;

        double pPrime = 1/(double)nPrime;

        parameterInput.get().setValue(idx, nPrime);

        return n*Math.log(1-pPrime) + Math.log(pPrime) - (nPrime*Math.log(1-p) + Math.log(p));
    }
}
