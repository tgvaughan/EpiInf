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

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class ScaleWithInt extends Operator {

    public Input<IntegerParameter> integerParameterInput = new Input<>(
            "integerParameter", "Scalar integer parameter to operate on.",
            Input.Validate.REQUIRED);

    public Input<Operator> operatorInput = new Input<>(
            "integerOperator", "Operator for integer parameter.",
            Input.Validate.REQUIRED);

    public Input<List<RealParameter>> realParameterInput = new Input<>(
            "realParameter", "Parameter to scale in proportion to the integer parameter.",
            new ArrayList<>());

    public Input<List<RealParameter>> realParameterInverseInput = new Input<>(
            "realParameterInverse", "Parameter to scale inversely to the integer parameter.",
            new ArrayList<>());

    Operator integerOperator;
    IntegerParameter integerParameter;

    @Override
    public void initAndValidate() {
        integerOperator = operatorInput.get();
        integerParameter = integerParameterInput.get();
    }

    @Override
    public double proposal() {

        double logHR;

        int n = integerParameter.getValue();

        if (n<=0)
            return Double.NEGATIVE_INFINITY;

        logHR = integerOperator.proposal();
        int m = integerParameter.getValue();

        if (m<=0)
            return Double.NEGATIVE_INFINITY;

        double f = m/(double)n;

        int upMinusDown = 0;
        for (RealParameter param : realParameterInput.get()) {
            for (int i=0; i<param.getDimension(); i++) {
                if(param.getValue(i) != 0.0) {
                    double newValue = param.getValue(i)*f;

                    if (newValue<param.getLower() || newValue>param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);

                    upMinusDown += 1;
                }
            }
        }

        for (RealParameter param : realParameterInverseInput.get()) {
            for (int i=0; i<param.getDimension(); i++) {
                if (param.getValue(i) != 0.0) {
                    double newValue = param.getValue(i)/f;

                    if (newValue<param.getLower() || newValue>param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);

                    upMinusDown -= 1;
                }
            }
        }

        logHR += Math.log(f)*upMinusDown;

        return logHR;
    }
}
