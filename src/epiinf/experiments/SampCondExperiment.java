package epiinf.experiments;

import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SABirthDeathModel;
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

        double tmin = 0.01;

        double dt = t0/nSteps;

        try (PrintStream ps = new PrintStream(args[0])) {

            ps.println("t1 densityEI densitySA");

            for (int i=0; i<nSteps; i++) {
                double t1 = dt*(i+1);

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

                RealParameter infectionRateParam = new RealParameter(String.valueOf(beta));
                RealParameter recoveryRateParam = new RealParameter(String.valueOf(gamma));
                RealParameter removalProbParam = new RealParameter("1.0");
                RealParameter rhoParam = new RealParameter("0.1");
                RealParameter psiParam = new RealParameter("0.0");
                RealParameter originParam = new RealParameter(String.valueOf(t0));

                BirthDeathModel model = new BirthDeathModel();
                model.initByName(
                        "infectionRate", infectionRateParam,
                        "recoveryRate", recoveryRateParam,
                        //"psiSamplingVariable", psiParam,
                        //"usePsiSamplingProportion", false,
                        "removalProb", removalProbParam,
                        //"rhoSamplingProb", rhoParam,
                        //"rhoSamplingTime", new RealParameter("0.0"),
                        //"rhoSamplingTimesBackward", true,
                        "origin", originParam);

                SMCTreeDensity density = new SMCTreeDensity();
                density.initByName(
                        "nParticles", 10000,
                        "useTauLeaping", false,
                        "model", model,
                        "finalSampleOffset", new RealParameter("0.0"),
                        "tree", tree);

                SABirthDeathModel densitySA = new SABirthDeathModel();
                densitySA.initByName(
                        "tree", tree,
                        "origin", originParam,
                        "birthRate", infectionRateParam,
                        "deathRate", recoveryRateParam,
                        "samplingRate", psiParam,
                        "rho", rhoParam,
                        "removalProbability", removalProbParam
                );

                ps.format("%g %g %g\n", t1,
                        density.calculateLogP(),
                        densitySA.calculateLogP());
            }

        } catch (FileNotFoundException ex) {

        }
    }
}
