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

package epiinf.models;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;

/**
 * General stochastic SIS model of an epidemic.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SISModel extends EpidemicModel {
    
    public Input<IntegerParameter> S0Input = new Input<>(
            "S0", "Initial size of susceptible population.", Validate.REQUIRED);

    @Override
    public EpidemicState getInitialState() {
        return new EpidemicState(S0Input.get().getValue(), 1, 0);
    }

    @Override
    protected double calculateInfectionPropensity(EpidemicState state) {
        return rateCache.get(state.modelIntervalIdx)[EpidemicEvent.INFECTION]*state.S*state.I;
    }

    @Override
    protected double calculateRecoveryPropensity(EpidemicState state) {
        return rateCache.get(state.modelIntervalIdx)[EpidemicEvent.RECOVERY]*state.I;
    }

    @Override
    public void incrementState(EpidemicState state, EpidemicEvent event) {
        switch(event.type) {
            case EpidemicEvent.INFECTION:
                state.S -= event.multiplicity;
                state.I += event.multiplicity;
                break;
            case EpidemicEvent.RECOVERY:
            case EpidemicEvent.RHO_SAMPLE:
            case EpidemicEvent.PSI_SAMPLE_REMOVE:
            case EpidemicEvent.OTHER_SAMPLE:
                state.I -= event.multiplicity;
                state.S += event.multiplicity;
                break;
            default:
                break;
        }
    }

    @Override
    public double getTau(double epsilon, EpidemicState state, double infectionProp, double recoveryProp) {
        double muS = -infectionProp + recoveryProp;
        double muI = infectionProp - recoveryProp;

        double sigmaS2 = infectionProp + recoveryProp;
        double sigmaI2 = infectionProp + recoveryProp;

        double epsS = Math.max(1.0, epsilon*state.S);
        double epsI = Math.max(1.0, epsilon*state.I);

        double tau = muS != 0.0 ? epsS/Math.abs(muS)
                : Double.POSITIVE_INFINITY;

        if (sigmaS2>0)
            tau = Math.min(tau, epsS*epsS/sigmaS2);

        if (muI>0)
            tau = Math.min(tau, epsI/Math.abs(muI));

        if (sigmaI2>0)
            tau = Math.min(tau, epsI*epsI/sigmaI2);

        return tau;
    }
}
