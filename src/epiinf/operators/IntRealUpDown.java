package epiinf.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * An up/down operator tailored for a product between a scalar integer parameter
 * and a vector real parameter.  It works by performing a random walk on the
 * integer and inducing a discrete random walk on the real.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IntRealUpDown extends Operator {

    public Input<IntegerParameter> intInput = new Input<>("int",
            "Integer parameter", Input.Validate.REQUIRED);

    public Input<RealParameter> realInput = new Input<>("real",
            "Real parameter", Input.Validate.REQUIRED);

    public Input<Integer> windowSizeInput = new Input<>("windowSize",
            "Window size for integer random walk.", Input.Validate.REQUIRED);

    private int windowSize;
    private RealParameter realParam;
    private IntegerParameter intParam;

    @Override
    public void initAndValidate() {
        intParam = intInput.get();
        realParam = realInput.get();
        windowSize = windowSizeInput.get();
    }

    @Override
    public double proposal() {

        final int newValue = intParam.getValue() + Randomizer.nextInt(2 * windowSize + 1) - windowSize;

        if (newValue < intParam.getLower() || newValue > intParam.getUpper())
            return Double.NEGATIVE_INFINITY;

        realParam.startEditing(this);
        realParam.scale(intParam.getArrayValue()/newValue);

        intParam.setValue(newValue);

        return 0;
    }

}
