package epiinf.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.evolution.tree.Tree;

import java.io.*;
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

    public int nCounts;
    public double[] ages;
    public int[] counts;

    public IncidenceData() { }

    @Override
    public void initAndValidate() {

        String string = null;
        if (valueInput.get() != null) {
            string = valueInput.get();
        } else if (fileNameInput.get() != null) {
            try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {
                if (fileHasHeaderInput.get())
                    reader.readLine();

                String nextLine;
                while ((nextLine = reader.readLine()) != null)
                    string += " " + nextLine;

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (string == null)
            throw new IllegalArgumentException("Must supply either fileName " +
                    "or value input.");
        String[] valueStrings = string.trim().split("\\s+");

        if (valueStrings.length % 2 != 0)
            throw new IllegalArgumentException("Error parsing");

        nCounts = valueStrings.length/2;

        ages = new double[nCounts];
        counts = new int[nCounts];

        for (int i=0; i<nCounts; i++) {
            ages[i] = Double.parseDouble(valueStrings[2*i]);
            counts[i] = Integer.parseInt(valueStrings[2*i+1]);
        }
    }
}
