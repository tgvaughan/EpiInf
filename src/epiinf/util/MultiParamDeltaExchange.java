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

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class MultiParamDeltaExchange extends Operator {

    public Input<List<RealParameter>> paramsInput = new Input<>(
            "realParameter",
            "Parameter to operate on.", new ArrayList<>());

    public Input<Double> relWindowSizeInput = new Input<>(
            "relWindowSize",
            "Size of window used to choose new parameter values, expressed " +
                    "as a fraction of the sum of all parameter values.", 0.1);

    List<RealParameter> params;
    double relWindowSize;

    @Override
    public void initAndValidate() {
        params = paramsInput.get();

        if (params.size()<2)
            throw new IllegalArgumentException(
                    "MultiParamDeltaExchange must be applied to two " +
                            "or more RealParameters.");

        relWindowSize = relWindowSizeInput.get();
    }

    @Override
    public double proposal() {

        int p1idx = Randomizer.nextInt(params.size());
        int p2idx = Randomizer.nextInt(params.size()-1);
        if (p2idx == p1idx)
            p2idx += 1;


        // Compute sum of all parameters:
        double sum = 0;
        for (int i=0; i<params.size(); i++)
            sum += params.get(i).getValue();

        double delta = sum*relWindowSize*Randomizer.nextDouble();
        double newP1 = params.get(p1idx).getValue() + delta;
        double newP2 = params.get(p2idx).getValue() - delta;

        if (newP1 > params.get(p1idx).getUpper()
                || newP2 < params.get(p2idx).getLower())
            return Double.NEGATIVE_INFINITY;

        params.get(p1idx).setValue(newP1);
        params.get(p2idx).setValue(newP2);

        return 0;
    }

}
