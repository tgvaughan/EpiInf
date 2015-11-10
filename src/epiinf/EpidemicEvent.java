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
    
    public static final EpidemicEvent Infection = new EpidemicEvent(-1, INFECTION, 1);
    public static final EpidemicEvent Recovery = new EpidemicEvent(-1, RECOVERY, 1);
    public static final EpidemicEvent RhoSample = new EpidemicEvent(-1, RHO_SAMPLE, 1);
    public static final EpidemicEvent PsiSampleRemove = new EpidemicEvent(-1, PSI_SAMPLE_REMOVE, 1);
    public static final EpidemicEvent PsiSampleNoRemove = new EpidemicEvent(-1, PSI_SAMPLE_NOREMOVE, 1);
    public static final EpidemicEvent OtherSample = new EpidemicEvent(-1, PSI_SAMPLE_REMOVE, 1);

    public static EpidemicEvent MultipleRhoSamples(int multiplicity) {
        return new EpidemicEvent(-1, RHO_SAMPLE, multiplicity);
    }

    public static EpidemicEvent MultipleOtherSamples(int multiplicity) {
        return new EpidemicEvent(-1, OTHER_SAMPLE, multiplicity);
    }
}
