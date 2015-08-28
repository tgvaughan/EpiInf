package epiinf.distribs;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import epiinf.models.SEISModel;

import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SEISTreeDensity extends TreeDistribution {

    public Input<SEISModel> modelInput = new Input<>(
            "model",
            "SEIS model under which to compute density.",
            Input.Validate.REQUIRED);

    TreeInterface tree;
    SEISModel model;

    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        tree = treeInput.get();
        model = modelInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        logP = 0.0;

        return logP;
    }
}

