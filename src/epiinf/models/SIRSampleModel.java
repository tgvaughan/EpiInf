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
public class SIRSampleModel extends EpidemicModel {
    
    public Input<IntegerParameter> S0Input = new Input<IntegerParameter>(
            "S0", "Initial size of susceptible population.", Validate.REQUIRED);
    
    public Input<RealParameter> infectionRateInput = new Input<RealParameter>(
            "infectionRate", "Infection rate.", Validate.REQUIRED);
    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
            "recoveryRate", "Recovery rate.", Validate.REQUIRED);
    
    public Input<RealParameter> samplingRateInput = new Input<RealParameter>(
            "samplingRate", "Rho sampling rate.", Validate.REQUIRED);

    @Override
    public EpidemicState getInitialState() {
        return new EpidemicState(S0Input.get().getValue(), 1, 0);
    }

    @Override
    public void calculatePropensities(EpidemicState state) {
        propensities.put(EpidemicEvent.EventType.INFECTION,
                infectionRateInput.get().getValue()*state.S*state.I);

        propensities.put(EpidemicEvent.EventType.RECOVERY,
                recoveryRateInput.get().getValue()*state.I);
        
        propensities.put(EpidemicEvent.EventType.SAMPLE,
                samplingRateInput.get().getValue()*state.I);

        totalPropensity = propensities.get(EpidemicEvent.EventType.INFECTION)
                + propensities.get(EpidemicEvent.EventType.RECOVERY);
    }

    @Override
    public void incrementState(EpidemicState state, EpidemicEvent.EventType type) {
        switch(type) {
            case INFECTION:
                state.S -= 1;
                state.I += 1;
                break;
            case RECOVERY:
                state.I -= 1;
                state.R += 1;
                break;
            case SAMPLE:
                state.I -= 1;
                state.R += 1;
            default:
                break;
        }
    }

    @Override
    public EpidemicEvent.EventType getCoalescenceEventType() {
        return EpidemicEvent.EventType.INFECTION;
    }

    @Override
    public EpidemicEvent.EventType getLeafEventType() {
        return EpidemicEvent.EventType.SAMPLE;
    }
    
}
