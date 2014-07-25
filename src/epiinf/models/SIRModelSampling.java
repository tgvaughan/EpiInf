/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package epiinf.models;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;

/**
 * Unstructured SIR model with an explicit (rho) sampling process.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SIRModelSampling extends EpidemicModel {
    
    public Input<IntegerParameter> S0Input = new Input<>(
            "S0", "Initial size of susceptible population.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<>(
            "infectionRate", "Infection rate.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<>(
            "recoveryRate", "Recovery rate.", Validate.REQUIRED);
    
    public Input<RealParameter> samplingProbInput = new Input<>(
            "samplingProb", "Rho sampling rate.", Validate.REQUIRED);

    @Override
    public EpidemicState getInitialState() {
        return new EpidemicState(S0Input.get().getValue(), 1, 0);
    }

    @Override
    public void calculatePropensities(EpidemicState state) {
        propensities.put(EpidemicEvent.Type.INFECTION,
                infectionRateInput.get().getValue()*state.S*state.I);

        propensities.put(EpidemicEvent.Type.RECOVERY,
                recoveryRateInput.get().getValue()
                        *(1.0-samplingProbInput.get().getValue())
                        *state.I);
        
        propensities.put(EpidemicEvent.Type.SAMPLE,
                recoveryRateInput.get().getValue()
                        *samplingProbInput.get().getValue()
                        *state.I);

        totalPropensity = propensities.get(EpidemicEvent.Type.INFECTION)
                + propensities.get(EpidemicEvent.Type.SAMPLE)
                + propensities.get(EpidemicEvent.Type.RECOVERY);
    }

    @Override
    public void incrementState(EpidemicState state, EpidemicEvent.Type type) {
        switch(type) {
            case INFECTION:
                state.S -= 1;
                state.I += 1;
                break;
            case RECOVERY:
            case SAMPLE:
            case MULTISAMPLE:
                state.I -= 1;
                state.R += 1;
                break;
            default:
                break;
        }
    }

    @Override
    public double getProbCoalescence(EpidemicState state, int lineages) {
        double N = state.I;
        
        if (N<2)
            return 0.0;
        else
            return 1.0/(N*(N-1));
    }

    @Override
    public double getProbNoCoalescence(EpidemicState state, int lineages) {
        double N = state.I;
        
        if (N<2)
            return 0.0;
        else
            return lineages*(lineages-1)/(N*(N-1));
    }
    
}
