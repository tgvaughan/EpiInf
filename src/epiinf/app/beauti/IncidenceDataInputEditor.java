package epiinf.app.beauti;

import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.inputeditor.InputEditor.ExpandOption;
import beastfx.app.util.FXUtils;
import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import epiinf.IncidenceData;
import javafx.application.Platform;
import javafx.embed.swing.SwingNode;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;

public class IncidenceDataInputEditor extends InputEditor.Base {

    IncidenceData incidenceData;
    JTextArea dataTextArea;
    JCheckBox timesAreAgesCheckbox;
    JComboBox<String> dateFormatBox;
    JTextField finalSampleOffsetTextField;

    JButton applyChangesButton;

    boolean loadingInProgress = false;

    public IncidenceDataInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return IncidenceData.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
        m_bAddButtons = addButtons;
        m_input = input;
        m_beastObject = beastObject;
        this.itemNr = itemNr;

        incidenceData = (IncidenceData) input.get();
        if (pane == null) {
        	pane = FXUtils.newHBox();
    		getChildren().add(pane);
        } else {
        	pane.getChildren().clear();
        }
    	addInputLabel();

        Box verticalBox = Box.createVerticalBox();
        verticalBox.setBorder(new EtchedBorder());

        Box horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(new JLabel(("Whitespace-delimited pairs of time and counts:")));
        horizontalBox.add(makeHorizontalFiller());
        verticalBox.add(horizontalBox);

        dataTextArea = new JTextArea(2, 30);
        dataTextArea.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                applyChangesButton.setEnabled(true);
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                applyChangesButton.setEnabled(true);
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                applyChangesButton.setEnabled(true);
            }
        });
        verticalBox.add(dataTextArea);

        applyChangesButton = new JButton("Accept changes");
        applyChangesButton.setEnabled(false);
        applyChangesButton.addActionListener(l -> {
            saveToModel();
        });
        JButton loadFromFileButton = new JButton("Load from file");
        loadFromFileButton.addActionListener(e -> {
           JFileChooser fileChooser = new JFileChooser();
           fileChooser.setDialogTitle("Select file containing incidence data...");
           int result = fileChooser.showDialog(null, "Load from file");
           if (result == JFileChooser.APPROVE_OPTION) {
               try (BufferedReader reader = new BufferedReader(new FileReader(fileChooser.getSelectedFile()))) {
                   StringBuilder sb = new StringBuilder();
                   reader.lines().forEachOrdered(l -> sb.append(l).append(" "));
                   dataTextArea.setText(sb.toString());

                   saveToModel();
               } catch (IOException fileNotFoundException) {
                   fileNotFoundException.printStackTrace();
               }
           }
        });
        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(applyChangesButton);
        horizontalBox.add(loadFromFileButton);
        horizontalBox.add(makeHorizontalFiller());
        verticalBox.add(horizontalBox);

        timesAreAgesCheckbox = new JCheckBox("Interpret times as ages before end of sampling period.");
        timesAreAgesCheckbox.addItemListener(l -> saveToModel());
        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(timesAreAgesCheckbox);
        horizontalBox.add(makeHorizontalFiller());
        verticalBox.add(horizontalBox);

        dateFormatBox = new JComboBox<>(new String[] {
                "Numeric",
                "yyyy-MM-dd",
                "dd-MM-yyyy",
                "MM-dd-yyyy",
                "yyyy/MM/dd",
                "dd/MM/yyyy",
                "MM/dd/yyyy"});
        dateFormatBox.setEditable(true);
        dateFormatBox.addActionListener(l -> saveToModel());
        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(new JLabel("Date format for times:"));
        horizontalBox.add(dateFormatBox);
        horizontalBox.add(makeHorizontalFiller());
        verticalBox.add(horizontalBox);

        horizontalBox = Box.createHorizontalBox();
        horizontalBox.add(new JLabel("Age of final incidence sample:"));
        finalSampleOffsetTextField = new JFormattedTextField(NumberFormat.getNumberInstance());
        finalSampleOffsetTextField.addActionListener(l -> saveToModel());
        horizontalBox.add(finalSampleOffsetTextField);
        horizontalBox.add(makeHorizontalFiller());
        verticalBox.add(horizontalBox);

        SwingNode n = new SwingNode();
        n.setContent(verticalBox);
        pane.getChildren().add(n);

        loadFromModel();
    }

    private void loadFromModel() {
        loadingInProgress = true;

        timesAreAgesCheckbox.setSelected(incidenceData.valuesAreAgesInput.get());
        dateFormatBox.setEnabled(!incidenceData.valuesAreAgesInput.get());
        finalSampleOffsetTextField.setEnabled(!incidenceData.valuesAreAgesInput.get());

        if (incidenceData.dateFormatInput.get() != null)
            dateFormatBox.setSelectedItem(incidenceData.dateFormatInput.get().trim());
        else
            dateFormatBox.setSelectedItem("Numeric");

        finalSampleOffsetTextField.setText(String.valueOf(incidenceData.finalSampleOffsetInput.get().getValue()));

        if (incidenceData.valueInput.get() != null)
            dataTextArea.setText(incidenceData.valueInput.get().trim());
        else
            dataTextArea.setText("");

        applyChangesButton.setEnabled(false);
        loadingInProgress = false;
    }

    private void saveToModel() {
        if (loadingInProgress)
            return;

        incidenceData.valuesAreAgesInput.setValue(timesAreAgesCheckbox.isSelected(), incidenceData);
        String dataString = dataTextArea.getText().trim();
        if (dataString.isEmpty())
            incidenceData.valueInput.setValue(null, incidenceData);
        else
            incidenceData.valueInput.setValue(dataString, incidenceData);

        RealParameter offsetParam = incidenceData.finalSampleOffsetInput.get();
        offsetParam.initByName("value", timesAreAgesCheckbox.isSelected()
                ? "0.0"
                : finalSampleOffsetTextField.getText());

        String dateFormat = (String) dateFormatBox.getSelectedItem();
        if (!timesAreAgesCheckbox.isSelected() && dateFormat != null && !dateFormat.equals("Numeric"))
            incidenceData.dateFormatInput.setValue(dateFormat, incidenceData);
        else
            incidenceData.dateFormatInput.setValue(null, incidenceData);

        refreshPanel();
        Platform.runLater(() -> 
        	init(m_input, m_beastObject, itemNr, ExpandOption.TRUE, m_bAddButtons)
        );
    }

    /**
     * @return a horizontal filler object.
     */
    public static Box.Filler makeHorizontalFiller() {
        return new Box.Filler(new Dimension(1,1),
                new Dimension(1,1),
                new Dimension(Integer.MAX_VALUE,1));
    }
}
