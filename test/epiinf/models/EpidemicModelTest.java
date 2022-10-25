package epiinf.models;

import beast.base.inference.parameter.RealParameter;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class EpidemicModelTest {

    @Test
    public void testBinarySearch() {

        RealParameter timeParameter = new RealParameter("20 50 100");

        RealParameter originParameter = new RealParameter("200.0");

        BirthDeathModel birthDeathModel = new BirthDeathModel();
        birthDeathModel.initByName(
                "origin", originParameter,
                "infectionRate", new RealParameter("2"),
                "recoveryRate", new RealParameter("1"),
                "removalProb", new RealParameter("0"));

        assertEquals(0, birthDeathModel.binarySearch(timeParameter, false, 0));
        assertEquals(0, birthDeathModel.binarySearch(timeParameter, false, 19));
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, false, 21));
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, false, 49));
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, false, 51));
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, false, 99));
        assertEquals(3, birthDeathModel.binarySearch(timeParameter, false, 101));

        assertEquals(3, birthDeathModel.binarySearch(timeParameter, true, 0));
        assertEquals(3, birthDeathModel.binarySearch(timeParameter, true, 99));
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, true, 101));
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, true, 149));
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, true, 151));
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, true, 179));
        assertEquals(0, birthDeathModel.binarySearch(timeParameter, true, 181));

        originParameter.setValue(101.0);
        assertEquals(3, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(99.0);
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(51.0);
        assertEquals(2, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(49.0);
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(21.0);
        assertEquals(1, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(19.0);
        assertEquals(0, birthDeathModel.binarySearch(timeParameter, true, 0));
        originParameter.setValue(0.0);
        assertEquals(0, birthDeathModel.binarySearch(timeParameter, true, 0));
    }
}

