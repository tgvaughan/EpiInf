/*
 * Copyright (C) 2016 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.core.BEASTObject;
import beast.core.Input;
import beast.evolution.tree.Tree;

import java.io.*;
import java.util.ArrayList;
import java.util.IllegalFormatException;
import java.util.List;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class IncidenceData extends BEASTObject {

    public Input<String> valueInput = new Input<>("value",
            "String containing pairs of whitespace-delimited " +
                    "age and incidence count pairs.");

    public Input<String> fileNameInput = new Input<>("fromFile",
            "Name of file containing pairs of whitespace-delimited " +
                    "age and count pairs.");

    public Input<Boolean> fileHasHeaderInput = new Input<>(
            "fileHasHeader",
            "If true, discards the first line of a provided input file.",
            false);

    public Input<Double> errorInput = new Input<>("error",
            "Maximum error in incidence times.", 0.0);

    private List<Double> ages = new ArrayList<>();

    public IncidenceData() { }

    @Override
    public void initAndValidate() {

        StringBuilder string = null;
        if (valueInput.get() != null) {
            string = new StringBuilder(valueInput.get());
        } else if (fileNameInput.get() != null) {
            string = new StringBuilder();

            boolean isHeader = fileHasHeaderInput.get();

            try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {
                String nextLine;
                while ((nextLine = reader.readLine()) != null) {
                    if (nextLine.startsWith("#"))
                        continue;

                    if (isHeader) {
                        isHeader = false;
                        continue;
                    }

                    string.append(" ").append(nextLine);
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (string == null)
            throw new IllegalArgumentException("Must supply either fileName " +
                    "or value input.");
        String[] valueStrings = string.toString().trim().split("\\s+");

        if (valueStrings.length % 2 != 0)
            throw new IllegalArgumentException("Error parsing");

        for (int i=0; i<valueStrings.length/2; i++) {
            double age = Double.parseDouble(valueStrings[2*i]);
            int count = Integer.parseInt(valueStrings[2*i + 1]);

            for (int j=0; j<count; j++)
                ages.add(age);
        }
    }

    public List<Double> getAges() {
        return ages;
    }

    public double getError() {
        return errorInput.get();
    }
}
