package epiinf;

import beast.base.core.Input;
import beast.base.core.Log;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;

public class TrajectoryFromFile extends EpidemicTrajectory {

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file containing trajectory.",
            Input.Validate.REQUIRED);

    public TrajectoryFromFile() { }

    @Override
    public void initAndValidate() {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileNameInput.get()))) {

            String header = reader.readLine();
            if (!header.equals("t S I R eventType eventMultiplicity"))
                throw new RuntimeException("Error parsing trajectory file header.  Is this a trajectory file?");

            String thisLine;
            while ((thisLine = reader.readLine()) != null) {
                String[] strVals = thisLine.split(" ");

                String eventType = strVals[4];
                if (eventType.equals("START")) {

                    EpidemicState state = new EpidemicState(
                            Double.valueOf(strVals[1]),
                            Double.valueOf(strVals[2]),
                            Double.valueOf(strVals[3]));
                    state.time = Double.valueOf(strVals[0]);

                    stateList.add(state);

                } else if (eventType.equals("END")) {

                    origin = Double.valueOf(strVals[0]);

                } else {

                    EpidemicEvent event = new EpidemicEvent();
                    event.setType(strVals[4])
                            .setTime(Double.valueOf(strVals[0]))
                            .setMultiplicity(Integer.valueOf(strVals[5]));

                    eventList.add(event);

                    EpidemicState state = new EpidemicState(
                            Double.valueOf(strVals[1]),
                            Double.valueOf(strVals[2]),
                            Double.valueOf(strVals[3]));
                    state.time = event.time;

                    stateList.add(state);
                }
            }

        } catch (java.io.IOException ex) {
            throw new RuntimeException("Error reading trajectory file '" + fileNameInput.get() + "'.");
        }
    }
}
