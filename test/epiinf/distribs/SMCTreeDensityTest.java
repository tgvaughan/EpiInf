/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import beast.util.TreeParser;
import epiinf.models.EpidemicModel;
import epiinf.models.SIRModel;
import epiinf.models.SISModel;

import static org.junit.Assert.assertTrue;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class SMCTreeDensityTest {

    @org.junit.Test
    public void testSISTreeDensityContemp() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "((t0:0.02075027867360646,t4:0.02075027867360646):3.0918239485" +
                        "84213,((t3:0.7877417253484267,t2:0.7877417253484267):" +
                        "0.6153647974911296,(t1:0.8814380219278748,t5:0.881438" +
                        "0219278748):0.5216685009116815):1.7094677044182633):0" +
                        ".8874257727421804;", false, false, true, 0);

        EpidemicModel model = new SISModel();
        model.initByName(
                "origin", new RealParameter("4.0"),
                "S0", new RealParameter("99"),
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.2"),
                "psiSamplingVariable", new RealParameter("0.0"),
                "rhoSamplingProb", new RealParameter("0.3"),
                "rhoSamplingTime", new RealParameter("4.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 500000);

        double logP = density.calculateLogP();
        double logPtrue = -5.85;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.01);
    }

    @org.junit.Test
    public void testSISTreeDensitySerial() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "(t19:0.5728982259951056,(t1:4.682548468426976,((((t13:0.1" +
                        "5140665947746434,t7:1.1107091305509993):1.8687945" +
                        "24893283,(((t3:1.2033210062102193,t0:1.5052311752" +
                        "862582):1.979155438506392,(t11:1.6751090237003017" +
                        ",t5:2.520636526974794):0.21768909221243238):0.383" +
                        "3431620211192,(t9:2.0664574891912832,t16:0.323283" +
                        "6660847038):0.5630070543827079):0.162464953866913" +
                        "4):0.4307319096416098,t14:2.0835649148745947):0.2" +
                        "091271957282972,(t2:4.036974668854562,((t17:0.620" +
                        "0967311580277,(t18:0.5326159975392679,(t12:1.7242" +
                        "66814761453,((t10:1.5764794559811675,t8:1.7898943" +
                        "65128709):0.5901428139973652,t15:0.81544003887419" +
                        "74):0.05238068412529273):0.007416803235933855):0." +
                        "008158146264380939):0.11241136563363074,(t6:2.254" +
                        "2045996975584,t4:2.809463362774489):0.56243066627" +
                        "73919):0.502771737926837):0.4410902670875849):0.0" +
                        "4714580222573783):0.06274869721383691):0.18596113" +
                        "70305343;", false, false, true, 0);

        EpidemicModel model = new SISModel();
        model.initByName(
                "origin", new RealParameter("4.96590947152"),
                "S0", new RealParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingVariable", new RealParameter("0.1"),
                "removalProb", new RealParameter("1.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 100000);

        double logP = density.calculateLogP();
        double logPtrue = -34.87;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.1);
    }

    @org.junit.Test
    public void testSISTreeDensitySerialLeap() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "(t19:0.5728982259951056,(t1:4.682548468426976,((((t13:0.1" +
                        "5140665947746434,t7:1.1107091305509993):1.8687945" +
                        "24893283,(((t3:1.2033210062102193,t0:1.5052311752" +
                        "862582):1.979155438506392,(t11:1.6751090237003017" +
                        ",t5:2.520636526974794):0.21768909221243238):0.383" +
                        "3431620211192,(t9:2.0664574891912832,t16:0.323283" +
                        "6660847038):0.5630070543827079):0.162464953866913" +
                        "4):0.4307319096416098,t14:2.0835649148745947):0.2" +
                        "091271957282972,(t2:4.036974668854562,((t17:0.620" +
                        "0967311580277,(t18:0.5326159975392679,(t12:1.7242" +
                        "66814761453,((t10:1.5764794559811675,t8:1.7898943" +
                        "65128709):0.5901428139973652,t15:0.81544003887419" +
                        "74):0.05238068412529273):0.007416803235933855):0." +
                        "008158146264380939):0.11241136563363074,(t6:2.254" +
                        "2045996975584,t4:2.809463362774489):0.56243066627" +
                        "73919):0.502771737926837):0.4410902670875849):0.0" +
                        "4714580222573783):0.06274869721383691):0.18596113" +
                        "70305343;", false, false, true, 0);

        EpidemicModel model = new SISModel();
        model.initByName(
                "S0", new RealParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "origin", new RealParameter("4.96590947152"),
                "psiSamplingVariable", new RealParameter("0.1"),
                "removalProb", new RealParameter("1.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 100000,
                "useTauLeaping", true);

        double logP = density.calculateLogP();
        double logPtrue = -34.87;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.5);
    }


    @org.junit.Test
    public void testSIRTreeDensityContemp() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "((t2:2.6540971886279987,(t6:0.5985184632462022,t3:0.598518" +
                        "4632462022):2.0555787253817965):1.3143407954121802" +
                        ",((((t11:1.7507486451524459,t18:1.7507486451524459" +
                        "):0.15463202333675685,(t7:0.005557722030708145,t19" +
                        ":0.005557722030708145):1.8998229464584946):0.04538" +
                        "965150569485,t17:1.9507703199948976):1.24654512383" +
                        "0562,((t4:1.3402939924855457,t0:1.3402939924855457" +
                        "):1.2267309476017108,(((t16:1.6297766592516973,(t8" +
                        ":1.604699304207743,((t14:0.7386842437142844,t9:0.7" +
                        "386842437142844):0.26795159631020393,t10:1.0066358" +
                        "400244884):0.5980634641832547):0.02507735504395425" +
                        "3):0.390258648420418,(t15:1.8472535458949886,(t20:" +
                        "0.9138121992788801,t22:0.9138121992788801):0.93344" +
                        "13466161084):0.17278176177712679):0.04086207977454" +
                        "6256,((t21:1.3594190980750436,(t13:1.0358015066819" +
                        "188,(t12:0.7248274453285872,t5:0.7248274453285872)" +
                        ":0.3109740613533316):0.32361759139312474):0.406905" +
                        "07105017915,t1:1.7663241691252227):0.2945732183214" +
                        "389):0.5061275526405948):0.6302905037382032):0.771" +
                        "1225402147193):0.03156201595982111;",
                false, false, true, 0);

        EpidemicModel model = new SIRModel();
        model.initByName(
                "origin", new RealParameter("4.0"),
                "S0", new RealParameter("199"),
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.2"),
                "psiSamplingVariable", new RealParameter("0.0"),
                "rhoSamplingProb", new RealParameter("0.3"),
                "rhoSamplingTime", new RealParameter("4.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 100000);

        double logP = density.calculateLogP();
        double logPtrue = -17.95;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.1);
    }

    @org.junit.Test
    public void testSIRTreeDensitySerial() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "((t10:1.190994645607128,((t5:0.7907784749295401,t0:2.137323" +
                        "582518187):1.4069968260134451,(((((t1:1.35747685412" +
                        "18845,t2:1.1326815198385916):1.353533897088965,t9:0" +
                        ".561536306309963):0.2839477233273353,t15:0.38749087" +
                        "865769694):0.08355380761559106,(t4:1.96826481583497" +
                        "83,t8:0.9918701157630321):0.18905207020723713):0.05" +
                        "1536457989408646,t7:1.3799674904273243):0.095042916" +
                        "80756127):0.13130794670059887):1.1041993932072272,(" +
                        "t13:1.796340301846171,(t12:1.7408704970223141,(((t1" +
                        "4:0.4222181780048908,t3:2.2744621113125665):0.36206" +
                        "299564941613,t11:1.124843163313686):0.1561269293807" +
                        "6744,t6:1.8854439746459617):0.4677724578417486):0.2" +
                        "1513248992333622):0.25084875675842433):0.1193998384" +
                        "1998814;" , false, false, true, 0);

        EpidemicModel model = new SIRModel();
        model.initByName(
                "origin", new RealParameter("4.89922758686"),
                "S0", new RealParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingVariable", new RealParameter("0.1"),
                "removalProb", new RealParameter("1.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 100000);

        double logP = density.calculateLogP();
        double logPtrue = -28.20;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.1);
    }

    @org.junit.Test
    public void testSIRTreeDensitySerialLeap() throws Exception {

        Randomizer.setSeed(42);

        TreeParser tree = new TreeParser(
                "((t10:1.190994645607128,((t5:0.7907784749295401,t0:2.137323" +
                        "582518187):1.4069968260134451,(((((t1:1.35747685412" +
                        "18845,t2:1.1326815198385916):1.353533897088965,t9:0" +
                        ".561536306309963):0.2839477233273353,t15:0.38749087" +
                        "865769694):0.08355380761559106,(t4:1.96826481583497" +
                        "83,t8:0.9918701157630321):0.18905207020723713):0.05" +
                        "1536457989408646,t7:1.3799674904273243):0.095042916" +
                        "80756127):0.13130794670059887):1.1041993932072272,(" +
                        "t13:1.796340301846171,(t12:1.7408704970223141,(((t1" +
                        "4:0.4222181780048908,t3:2.2744621113125665):0.36206" +
                        "299564941613,t11:1.124843163313686):0.1561269293807" +
                        "6744,t6:1.8854439746459617):0.4677724578417486):0.2" +
                        "1513248992333622):0.25084875675842433):0.1193998384" +
                        "1998814;" , false, false, true, 0);

        EpidemicModel model = new SIRModel();
        model.initByName(
                "origin", new RealParameter("4.89922758686"),
                "S0", new RealParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingVariable", new RealParameter("0.1"),
                "removalProb", new RealParameter("1.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "model", model,
                "finalSampleOffset", new RealParameter("0.0"),
                "nParticles", 100000,
                "useTauLeaping", true);

        double logP = density.calculateLogP();
        double logPtrue = -28.20;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.1);
    }

}
