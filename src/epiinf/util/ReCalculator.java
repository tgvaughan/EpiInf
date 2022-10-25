package epiinf.util;

import beast.base.core.BEASTObject;
import beast.base.core.Function;
import beast.base.core.Input;

import java.util.ArrayList;
import java.util.List;

public class ReCalculator extends BEASTObject implements Function {

    public Input<Function> infectionRateInput = new Input<>("infectionRate", "Infection rate parameter", Input.Validate.REQUIRED);
    public Input<Function> infectionRateChangeTimesInput = new Input<>("infectionRateChangeTimes", "Infection rate change times");
    public Input<Function> deathRateInput = new Input<>("deathRate", "Death rate parameter", Input.Validate.REQUIRED);
    public Input<Function> deathRateChangeTimesInput = new Input<>("deathRateChangeTimes", "Infection rate change times");
    public Input<Function> samplingPropInput = new Input<>("samplingProp", "Sampling proportion parameter", Input.Validate.REQUIRED);
    public Input<Function> samplingProbChangeTimesInput = new Input<>("samplingPropChangeTimes", "Infection rate change times");

    Function infectionRate, deathRate, samplingProp;
    Function infectionRateChangeTimes, deathRateChangeTimes, samplingPropChangeTimes;

    List<Double> Re, ReChangeTimes;

    @Override
    public void initAndValidate() {
        infectionRate = infectionRateInput.get();
        deathRate = deathRateInput.get();
        samplingProp = samplingPropInput.get();

        infectionRateChangeTimes = infectionRateChangeTimesInput.get();
        deathRateChangeTimes = deathRateChangeTimesInput.get();
        samplingPropChangeTimes = samplingProbChangeTimesInput.get();

        Re = new ArrayList<Double>();
        ReChangeTimes = new ArrayList<Double>();

    }

    void update() {
        int infIdx=0, deathIdx=0, sampIdx=0;

        Re.clear();
        ReChangeTimes.clear();

        double nextChangeTime = Double.POSITIVE_INFINITY;

        do {
            Re.add(infectionRate.getArrayValue(infIdx)/deathRate.getArrayValue(deathIdx)*(1-samplingProp.getArrayValue(sampIdx)));

            if (infectionRateChangeTimes != null && infIdx < infectionRateChangeTimes.getDimension())
                nextChangeTime = Math.min(nextChangeTime, infectionRateChangeTimes.getArrayValue(deathIdx));
            if (deathRateChangeTimes != null && deathIdx < deathRateChangeTimes.getDimension())
                nextChangeTime = Math.min(nextChangeTime, deathRateChangeTimes.getArrayValue(deathIdx));
            if (samplingPropChangeTimes != null && sampIdx < samplingPropChangeTimes.getDimension())
                nextChangeTime = Math.min(nextChangeTime, samplingPropChangeTimes.getArrayValue(sampIdx));

            if (infectionRateChangeTimes != null && infIdx < infectionRateChangeTimes.getDimension()
                    && nextChangeTime == infectionRateChangeTimes.getArrayValue(infIdx))
                infIdx += 1;

            if (deathRateChangeTimes != null && deathIdx < deathRateChangeTimes.getDimension()
                    && nextChangeTime == deathRateChangeTimes.getArrayValue(deathIdx))
                deathIdx += 1;

            if (samplingPropChangeTimes != null && sampIdx < samplingPropChangeTimes.getDimension()
                    && nextChangeTime == samplingPropChangeTimes.getArrayValue(sampIdx))
                sampIdx += 1;

            if (nextChangeTime < Double.POSITIVE_INFINITY)
                ReChangeTimes.add(nextChangeTime);

        } while(nextChangeTime<Double.POSITIVE_INFINITY);
    }

    @Override
    public int getDimension() {
        update();

        return Re.size();
    }

    @Override
    public double getArrayValue(int dim) {
        update();

        return Re.get(dim);
    }
}
