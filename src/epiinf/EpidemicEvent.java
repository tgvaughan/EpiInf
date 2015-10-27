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
    public enum Type { INFECTION, RECOVERY, RHO_SAMPLE, PSI_SAMPLE };

    public Type type;
    public int multiplicity;

    public EpidemicEvent() {
        multiplicity = 1;
    }
    
    public EpidemicEvent(double time, Type type, int multiplicity) {
        this.time = time;
        this.type = type;
        this.multiplicity = multiplicity;
    }
    
    /**
     * Get record describing this event for use in constructing
     * R-compatible text files.
     * 
     * @return record
     */
    public String getRecord() {
        return type + " " + multiplicity;
    }

    @Override
    public String toString() {
        return "t: " + time + " type: " + type + " mult: " + multiplicity;
    }
    
    public static EpidemicEvent Infection = new EpidemicEvent(-1, Type.INFECTION, 1);
    public static EpidemicEvent Recovery = new EpidemicEvent(-1, Type.RECOVERY, 1);
    public static EpidemicEvent RhoSample = new EpidemicEvent(-1, Type.RHO_SAMPLE, 1);
    public static EpidemicEvent PsiSample = new EpidemicEvent(-1, Type.PSI_SAMPLE, 1);

    public static EpidemicEvent MultipleRhoSamples(int multiplicity) {
        return new EpidemicEvent(-1, Type.RHO_SAMPLE, multiplicity);
    }
}
