package epiinf.experiments;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import epiinf.distribs.SMCTreeDensity;
import epiinf.models.BirthDeathModel;

import java.io.FileNotFoundException;
import java.io.PrintStream;

public class SampCondExperiment {

    public static void main(String[] args) {

        if (args.length == 0)
            throw new IllegalArgumentException("File name command line argument needed.");

        int nSteps = 100;
        double t0 = 2.0;
        double beta = 2.0;
        double gamma = 0.0;

        double tmin = 0.1;

        double dt = (t0 - tmin)/(nSteps-1);

        try (PrintStream ps = new PrintStream(args[0])) {

            ps.println("t1 density");

            for (int i=0; i<nSteps; i++) {
                double t1 = tmin + dt*i;

                Node root = new Node();
                root.setHeight(t1);
                root.setNr(2);

                Node child1 = new Node("t1");
                child1.setHeight(0.0);
                child1.setNr(0);
                root.addChild(child1);

                Node child2 = new Node("t2");
                child2.setHeight(0.0);
                child2.setNr(1);
                root.addChild(child2);

                Tree tree = new Tree(root);
                tree.initAndValidate();

                BirthDeathModel model = new BirthDeathModel();
                model.initByName(
                        "infectionRate", new RealParameter(String.valueOf(beta)),
                        "recoveryRate", new RealParameter(String.valueOf(gamma)),
                        "removalProb", new RealParameter("1.0"),
                        "rhoSamplingProb", new RealParameter("0.0"),
                        "origin", new RealParameter(String.valueOf(t0)));

                SMCTreeDensity density = new SMCTreeDensity();
                density.initByName(
                        "nParticles", 10000,
                        "useTauLeaping", false,
                        "model", model,
                        "finalSampleOffset", new RealParameter("0.0"),
                        "tree", tree);

                ps.format("%g %g\n", t1, Math.exp(density.calculateLogP()));
            }

        } catch (FileNotFoundException ex) {

        }
    }
}
