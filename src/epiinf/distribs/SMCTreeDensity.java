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

package epiinf.distribs;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.Binomial;
import beast.util.Randomizer;
import beast.util.TreeParser;
import com.google.common.collect.Lists;
import epiinf.EpidemicEvent;
import epiinf.EpidemicState;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import epiinf.models.EpidemicModel;
import epiinf.models.SISModel;
import java.io.PrintStream;
import java.util.List;
import java.util.Random;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
@Description("Use SMC to estimate density of tree conditional on model "
        + "parameters.")
public class SMCTreeDensity extends Distribution {

    public Input<EpidemicModel> modelInput = new Input<>(
            "model", "Epidemic model.", Validate.REQUIRED);
    
    public Input<TreeEventList> treeEventListInput = new Input<>(
            "treeEventList", "Tree event list.", Validate.REQUIRED);
    
    public Input<Integer> nParticlesInput = new Input<>(
            "nParticles", "Number of particles to use in SMC calculation.",
            Validate.REQUIRED);

    EpidemicModel model;
    TreeEventList eventList;
    int nParticles;
    
    // DEBUG
//    PrintStream debugOut;

    public SMCTreeDensity() { }

    @Override
    public void initAndValidate() throws Exception {
        model = modelInput.get();
        eventList = treeEventListInput.get();
        nParticles = nParticlesInput.get();
    }

    @Override
    public double calculateLogP() throws Exception {
        
        // DEBUG
//        debugOut = new PrintStream("SMCdebug.json");
//        debugOut.println("{");
        
        logP = 0.0;
        
        List<Double> particleWeights = Lists.newArrayList();
        List<EpidemicState> particleStates = Lists.newArrayList();
        List<EpidemicState> particleStatesNew = Lists.newArrayList();
        
        // Initialize particles
        for (int p=0; p<nParticles; p++)
            particleStates.add(model.getInitialState());
        
        double t = 0.0;
        int k = 1;
//        int interval = 0;
        for (TreeEvent treeEvent : eventList.getEventList()) {
            
            // DEBUG
//            if (interval>0)
//                debugOut.println(",");
//            debugOut.println("\"interval" + interval + "\": {");
            
            // Update particles
            particleWeights.clear();
            double sumOfWeights = 0.0;
            for (int p=0; p<nParticles; p++) {
                
                // DEBUG
//                if (p>0)
//                    debugOut.println(", ");
//                debugOut.println("\"p" + p + "\": {");
                
                double newWeight = updateParticle(particleStates.get(p), t, k, treeEvent);
                
                particleWeights.add(newWeight);                
                sumOfWeights += newWeight;
                
                // DEBUG
//                debugOut.print("\n}");
            }
            
            // DEBUG
//            debugOut.print("}");
            
            // Update marginal likelihood estimate
            logP += Math.log(sumOfWeights/nParticles);
            
            if (!(sumOfWeights>0.0))
                return Double.NEGATIVE_INFINITY;

            // Sample particle with replacement
            particleStatesNew.clear();
            for (int p=0; p<nParticles; p++) {
                double u = Randomizer.nextDouble()*sumOfWeights;
                
                int pChoice;
                for (pChoice = 0; pChoice<nParticles; pChoice++) {
                    u -= particleWeights.get(pChoice);
                    if (u<0.0)
                        break;
                }
                
                if (pChoice == nParticles)
                    System.err.println("sumOfWeights: " + sumOfWeights);
                
                particleStatesNew.add(particleStates.get(pChoice).copy());
            }
            
            // Switch particleStates and particleStatesNew
            List<EpidemicState> temp = particleStates;
            particleStates = particleStatesNew;
            particleStatesNew = temp;
            
            // Update lineage counter
            if (treeEvent.type == TreeEvent.Type.COALESCENCE)
                k += 1;
            else
                k -= treeEvent.multiplicity;
            
            // Update start interval time
            t = treeEvent.time;
            
//            interval += 1;
        }
        
        // DEBUG
//        debugOut.println("}");
        
        return logP;
    }

    /**
     * Updates weight and state of particle.
     * 
     * @param particleState
     * @param startTime
     * @param lineages
     * @param finalTreeEvent
     * @return conditional prob of tree interval under trajectory
     */
    private double updateParticle(EpidemicState particleState,
            double startTime, int lineages, TreeEvent finalTreeEvent) {
        double conditionalP = 1.0;
        
        double t = startTime;
        
        // DEBUG
//        List<Double> tList = Lists.newArrayList();
//        List<Double> nList = Lists.newArrayList();
//        tList.add(t);
//        nList.add(particleState.I);
        
        while (true) {
            model.calculatePropensities(particleState);

            if (model.getTotalPropensity() > 0.0)
                t += Randomizer.nextExponential(model.getTotalPropensity());
            else
                t = Double.POSITIVE_INFINITY;
            
            if (t>finalTreeEvent.time)
                break;
            
            // Check for missed rho sampling event
            if (t > model.getNextRhoSamplingTime(t))
                return 0.0;
            
            double u = model.getTotalPropensity()*Randomizer.nextDouble();
            
            EpidemicEvent event = new EpidemicEvent();
            for (EpidemicEvent.Type type : model.getPropensities().keySet()) {
                u -= model.getPropensities().get(type);
                
                if (u<0) {
                    event.type = type;
                    break;
                }
            }
            
            if (event.type == EpidemicEvent.Type.RECOVERY) {
                conditionalP *= 1.0 - lineages/particleState.I; // prob not on sampled lineage

		// Prob that given recovery is not on sampled lineage, sampling
		// did not occur.
                if (model.psiSamplingProbInput.get() != null)
                    conditionalP *= 1.0 - model.psiSamplingProbInput.get().getValue();
            }

            model.incrementState(particleState, event);

            // Early exit if invalid state:
            if (particleState.I<lineages)
                return 0.0;
            
            // Increment conditional prob
            if (event.type == EpidemicEvent.Type.INFECTION)
                conditionalP *= model.getProbNoCoalescence(particleState, lineages);

            
            // DEBUG
//            tList.add(t);
//            nList.add(particleState.I);
            
        }
        
        // Include probability of tree event
        if (finalTreeEvent.type == TreeEvent.Type.COALESCENCE) {
            model.incrementState(particleState, EpidemicEvent.Infection);
            model.calculatePropensities(particleState);
            conditionalP *= model.getProbCoalescence(particleState, lineages+1)
                    *model.getPropensities().get(EpidemicEvent.Type.INFECTION);
        } else {
            
            // If the model contains an explicit sampling process, evaluate the
            // probability of the sampling event on the tree
            if (!model.rhoSamplingProbInput.get().isEmpty()
                    || model.psiSamplingProbInput.get() != null) {
            
                double sampleProb = 0.0;
                
                // If model contains a rho sampling event at this time, calculate the probability
                // of sampling the number of samples in finalTreeEvent given the current
                // state.
                for (int i=0; i<model.rhoSamplingProbInput.get().size(); i++) {
                    double rhoProb = model.rhoSamplingProbInput.get().get(i).getValue();
                    double rhoTime = model.rhoSamplingTimeInput.get().get(i).getValue();
                    
                    if (Math.abs(rhoTime - finalTreeEvent.time)<model.getTolerance()) {
                        int I = (int)Math.round(particleState.I);
                        int k = finalTreeEvent.multiplicity;
                        sampleProb += Binomial.choose(I, k)
                                *Math.pow(rhoProb, k)*Math.pow(1.0-rhoProb, I-k);
                    }
                }
                
                // If the model contains a non-zero psi sampling rate, calculate the
                // probability of a sampled recovery occuring at the time of finalTreeEvent
                
                if (finalTreeEvent.multiplicity==1
                        && model.psiSamplingProbInput.get() != null) {
                    model.calculatePropensities(particleState);
                    sampleProb += model.getPropensities().get(EpidemicEvent.Type.RECOVERY)
                            *model.psiSamplingProbInput.get().getValue();
                }
                
                conditionalP *= sampleProb;//*Math.exp(GammaFunction.lnGamma(1+finalTreeEvent.multiplicity));

            }
            
            model.incrementState(particleState,
                    EpidemicEvent.MultipleSamples(finalTreeEvent.multiplicity));

        }
        
        // DEBUG
//        tList.add(finalTreeEvent.time);
//        nList.add(particleState.I);
//        debugOut.print("\"t\": [");
//        for (int s=0; s<tList.size(); s++) {
//            if (s>0)
//                debugOut.print(",");
//            debugOut.print(tList.get(s));
//        }
//        debugOut.print("], \"n\": [");
//        for (int s=0; s<nList.size(); s++) {
//            if (s>0)
//                debugOut.print(",");
//            debugOut.print(nList.get(s));
//        }
//        debugOut.print("]");

        if (!particleState.isValid())
            return 0.0;
        else
            return conditionalP;
    }
    
    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
    
    /**
     * Main method for testing/debugging.
     * 
     * @param args 
     * @throws java.lang.Exception 
     */
    public static void main (String [] args) throws Exception {
        
        //Randomizer.setSeed(42);
        
        SISModel model = new SISModel();
//        model.initByName(
//                "S0", new IntegerParameter("99"),
//                "infectionRate", new RealParameter("0.01"),
//                "recoveryRate", new RealParameter("0.2"));
        
//        TrajectorySimulator trajSim = new TrajectorySimulator();
//        trajSim.initByName(
//                "model", model,
//                "maxDuration", 50.0,
//                "samplingTime", 12.0,
//                "sampleSize", 1000,
//                "allowOversampling", true,
//                "fileName", "truth.txt");
//        
        //Tree tree = new Tree();
//        TreeParser tree = new TreeParser(
//                Files.toString(new File("tims_tree.newick"),
//                        Charset.defaultCharset()));
        TreeParser tree = new TreeParser("((((((((54:0.760457025379198,(71:0.3437653858620582,45:0.3437653858620582)88:0.41669163951713983)97:2.474504945642595,((33:0.22603166792133855,1:0.22603166792133855)85:0.30931274992853375,48:0.5353444178498723)92:2.6996175531719206)121:0.20881179211241552,(62:0.02595425535515794,46:0.02595425535515794)78:3.4178195077790505)126:0.4721680834751627,((77:1.6361436522189265,(14:0.9873792093350531,22:0.9873792093350531)99:0.6487644428838735)108:0.3184342009557124,68:1.954577853174639)110:1.9613639934347322)130:2.6167030830107665,(4:1.2581270629966461,64:1.2581270629966461)103:5.2745178666234915)149:1.0127382692963423,(((12:0.041072598449853004,20:0.041072598449853004)79:2.8358789249741783,(32:0.48795777428765064,60:0.48795777428765064)91:2.3889937491363806)117:3.7409358636461736,((((55:1.4571583187659982,(63:1.302517704430663,34:1.302517704430663)104:0.15464061433533516)106:2.296218949145983,51:3.753377267911981)128:0.8065529497682551,(29:2.8831257023612995,43:2.8831257023612995)118:1.6768045153189366)135:1.405501594523245,((((19:2.207130466714924,44:2.207130466714924)113:1.75510934483348,10:3.962239811548404)131:0.8809458098515037,36:4.843185621399908)137:0.7862471727125779,42:5.6294327941124855)142:0.3359990180909955)145:0.6524555748667238)150:0.927495811846275)152:0.11318080604516556,(((((7:3.271502511116257,38:3.271502511116257)124:1.1332999977293925,(15:3.2397158923520575,31:3.2397158923520575)122:1.1650866164935918)134:1.2596142049326975,(58:5.311608185574999,(69:0.36602016135929283,6:0.36602016135929283)89:4.945588024215706)141:0.3528085282033482)143:0.03915719171337617,(((49:0.3167594954859947,59:0.3167594954859947)87:2.237346290073294,(((40:1.0756924007215378,((57:0.07076529412625554,2:0.07076529412625554)82:0.5690125241767827,(74:0.48388643224580186,17:0.48388643224580186)90:0.1558913860572364)93:0.43591458241849956)101:0.2561193641195256,(47:0.6404674774833179,(5:0.06345364386015184,65:0.06345364386015184)81:0.577013833623166)95:0.6913442873577456)105:0.15610896879505098,30:1.4879207336361144)107:1.0661850519231741)115:0.716140695075616,41:3.2702464806349045)123:2.4333274248568184)144:1.3167141933052244,((0:5.245522850006264,(37:3.151882105988957,(9:1.9589888304643797,(75:1.1489114505705604,76:1.1489114505705604)102:0.8100773798938192)111:1.1928932755245771)120:2.093640744017307)140:0.9661343340944768,((((((52:3.295386949497356,(72:1.6745218383528186,13:1.6745218383528186)109:1.6208651111445374)125:0.5951066972135521,(70:2.9204238091372456,(28:2.1979703746349912,61:2.1979703746349912)112:0.7224534345022544)119:0.9700698375736625)129:0.27526927078601204,24:4.16576291749692)132:0.6267752062569878,73:4.792538123753908)136:1.3461730402299255,((66:0.9600862830787236,39:0.9600862830787236)98:1.7138455443812006,((26:0.23744322469280377,18:0.23744322469280377)86:0.4207826951031848,3:0.6582259197959885)96:2.0157059076639356)116:3.4647793365239092)146:0.049851489369084234,((((27:0.04314073518095718,25:0.04314073518095718)80:2.224023195570995,8:2.267163930751952)114:1.3453438646189113,56:3.6125077953708633)127:1.4799551124012975,(35:0.6399728736472543,21:0.6399728736472543)94:4.4524900341249065)139:1.0960997455807568)147:0.023094530747822972)148:0.8086309146962067)151:0.6382759061646981)153:3.1840316089012273,(((67:1.0107319187570187,23:1.0107319187570187)100:3.161465828804345,11:4.172197747561364)133:0.7779517346178721,(50:0.20453084663303578,(53:0.09996733089231746,16:0.09996733089231746)83:0.10456351574071832)84:4.7456186355462)138:5.892446131683637)154:0.0");
        RealParameter treeOrigin = new RealParameter("12.0");
        
//        TransmissionTreeSimulator treeSim = new TransmissionTreeSimulator();
//        treeSim.initByName(
//                "tree", tree,
//                "treeOrigin", treeOrigin,
//                "epidemicTrajectory", trajSim,
//                "model", model,
//                "fileName", "truth.newick");
//        treeSim.initStateNodes();
        
        TreeEventList treeEventList = new TreeEventList();
        treeEventList.initByName(
                "tree", tree,
                "treeOrigin", treeOrigin); 
        
        try (PrintStream ps = new PrintStream("tims_tree.txt")) {
            treeEventList.writeExpoTreeFile(ps);
        }
        
        SMCTreeDensity treeDensity = new SMCTreeDensity();
        
        try (PrintStream ps = new PrintStream("logLik.txt")) {
            ps.println("beta logP");
            for (double beta=0.005; beta<0.015; beta += 0.0005) {
                //double beta = 0.01*Math.pow(1.3, bidx);
                //double beta = 0.01;
                
                model.initByName(
                        "S0", new IntegerParameter("99"),
                        "infectionRate", new RealParameter(String.valueOf(beta)),
                        "recoveryRate", new RealParameter("0.2"),
                        "rhoSamplingProb", new RealParameter("1.0"),
                        "rhoSamplingTime", new RealParameter("12.0"));
                
                treeDensity.initByName(
                        "treeEventList", treeEventList,
                        "model", model,
                        "nParticles", 4000);
                
                double logP = treeDensity.calculateLogP();
                
                System.out.println("beta: " + beta
                        + " logP: " + logP);
                ps.println(beta + " " + logP);
            }
        }
    }
}
