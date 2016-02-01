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
 * A state of an epidemic trajectory.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpidemicState {
    public double S, I, R;
    public double time;
    public int modelIntervalIdx, treeIntervalIdx;

    public EpidemicState() { }
    
    public EpidemicState(double S, double I, double R) {
        this.S = S;
        this.I = I;
        this.R = R;
        this.time = 0;
        this.modelIntervalIdx = 0;
        this.treeIntervalIdx = 0;
    }

    /**
     * Test whether state is valid or not.
     * 
     * @return true if state is valid.
     */
    public boolean isValid() {
        return (this.S>=0 && this.I>=0 && this.R>=0);
    }
    
    public EpidemicState copy() {
        EpidemicState stateCopy = new EpidemicState();
        stateCopy.S = S;
        stateCopy.I = I;
        stateCopy.R = R;
        stateCopy.time = time;
        stateCopy.modelIntervalIdx = modelIntervalIdx;
        stateCopy.treeIntervalIdx = treeIntervalIdx;

        return stateCopy;
    }

    public void assignFrom(EpidemicState otherState) {
        S = otherState.S;
        I = otherState.I;
        R = otherState.R;
        time = otherState.time;
        modelIntervalIdx = otherState.modelIntervalIdx;
        treeIntervalIdx = otherState.treeIntervalIdx;
    }

    @Override
    public String toString() {
        return "S: " + S + ", I: " + I + ", R: " + R;
    }
    
    /**
     * Retrieve header suitable for dumping rows of states to an R-compatible
     * file.
     * 
     * @return R input file header
     */
    public static String getHeader() {
        return "S I R";
    }
    
    /**
     * Obtain this state as a record in an R-compatible text file.
     * 
     * @return R input file record.
     */
    public String getRecord() {
        return (long)S + " " + (long)I + " " + (long)R;
    }
}
