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

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import beast.util.TreeParser;
import epiinf.OriginFromTrajectory;
import epiinf.models.EpidemicModel;
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
                "S0", new IntegerParameter("99"),
                "infectionRate", new RealParameter("0.01"),
                "recoveryRate", new RealParameter("0.2"),
                "psiSamplingRate", new RealParameter("0.0"),
                "rhoSamplingProb", new RealParameter("0.3"),
                "rhoSamplingTime", new RealParameter("4.0"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "treeOrigin", new RealParameter("4.0"),
                "model", model,
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
                "S0", new IntegerParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingRate", new RealParameter("0.1"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "treeOrigin", new RealParameter("4.96590947152"),
                "model", model,
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
                "S0", new IntegerParameter("99"),
                "infectionRate", new RealParameter("0.02"),
                "recoveryRate", new RealParameter("0.1"),
                "psiSamplingRate", new RealParameter("0.1"));

        SMCTreeDensity density = new SMCTreeDensity();
        density.initByName(
                "tree", tree,
                "treeOrigin", new RealParameter("4.96590947152"),
                "model", model,
                "nParticles", 100000,
                "nTauLeaps", 10,
                "alpha", 0.0);

        double logP = density.calculateLogP();
        double logPtrue = -34.87;

        System.out.println("Truth: " + logPtrue);
        System.out.println("Estimate: " + logP);

        assertTrue(Math.abs(logP-logPtrue)<0.5);
    }
}
