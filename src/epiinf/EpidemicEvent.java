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

package epiinf;

import java.util.Arrays;

/**
 * Class representing events during an epidemic.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpidemicEvent extends Event {
//    public enum Type { INFECTION, RECOVERY, RHO_SAMPLE, PSI_SAMPLE_REMOVE, OTHER_SAMPLE};
    public static final int INFECTION = 0;
    public static final int RECOVERY = 1;
    public static final int RHO_SAMPLE = 2;
    public static final int PSI_SAMPLE_REMOVE = 3;
    public static final int PSI_SAMPLE_NOREMOVE = 4;
    public static final int OTHER_SAMPLE= 5;
    public static final int nTypes = 6;
    public static final String[] typeNames = {
            "INFECTION",
            "RECOVERY",
            "RHO_SAMPLE",
            "PSI_SAMPLE_REMOVE",
            "PSI_SAMPLE_NOREMOVE",
            "OTHER_SAMPLE"
    };

    public int type;
    public int multiplicity;

    public EpidemicEvent() {
        multiplicity = 1;
    }
    
    public EpidemicEvent(double time, int type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }

    /**
     * @return true iff this is a sampling event
     */
    public boolean isSample() {
        return type >= 2;
    }

    /**
     * Get record describing this event for use in constructing
     * R-compatible text files.
     * 
     * @return record
     */
    public String getRecord() {
        return typeNames[type] + " " + multiplicity;
    }

    @Override
    public String toString() {
        return "t: " + time + " type: " + typeNames[type] + " mult: " + multiplicity;
    }

    /**
     * @return string representation of event type
     */
    public String getTypeName() {
        return typeNames[type];
    }

    /**
     * Set type of event.
     *
     * @param type new type
     * @return object reference for method chaining.
     */
    public EpidemicEvent setType(int type) {
        this.type = type;
        return this;
    }

    /**
     * Set type of event corresponding to type name.
     *
     * @param typeName name of type
     * @return object reference for method chaining.
     */
    public EpidemicEvent setType(String typeName) {
        type = Arrays.asList(typeNames).indexOf(typeName);
        return this;
    }

    /**
     * Set time of event.
     *
     * @param time new time
     * @return object reference for method chaining.
     */
    public EpidemicEvent setTime(double time) {
        this.time = time;
        return this;
    }

    /**
     * Set multiplicity of event.
     *
     * @param multiplicity new multiplicity
     * @return object reference for method chaining.
     */
    public EpidemicEvent setMultiplicity(int multiplicity) {
        this.multiplicity = multiplicity;
        return this;
    }

    public static final EpidemicEvent Infection = new EpidemicEvent(-1, INFECTION, 1);
    public static final EpidemicEvent Recovery = new EpidemicEvent(-1, RECOVERY, 1);
    public static final EpidemicEvent RhoSample = new EpidemicEvent(-1, RHO_SAMPLE, 1);
    public static final EpidemicEvent PsiSampleRemove = new EpidemicEvent(-1, PSI_SAMPLE_REMOVE, 1);
    public static final EpidemicEvent PsiSampleNoRemove = new EpidemicEvent(-1, PSI_SAMPLE_NOREMOVE, 1);
    public static final EpidemicEvent OtherSampleRemove = new EpidemicEvent(-1, PSI_SAMPLE_REMOVE, 1);
    public static final EpidemicEvent OtherSampleNoRemove = new EpidemicEvent(-1, PSI_SAMPLE_NOREMOVE, 1);

    public static EpidemicEvent MultipleInfections(int multiplicity) {
        return new EpidemicEvent(-1, INFECTION, multiplicity);
    }

    public static EpidemicEvent MultipleRecoveries(int multiplicity) {
        return new EpidemicEvent(-1, RECOVERY, multiplicity);
    }

    public static EpidemicEvent MultiplePsiSampleRemove(int multiplicity) {
        return new EpidemicEvent(-1, RECOVERY, multiplicity);
    }


    public static EpidemicEvent MultipleRhoSamples(int multiplicity) {
        return new EpidemicEvent(-1, RHO_SAMPLE, multiplicity);
    }
}
