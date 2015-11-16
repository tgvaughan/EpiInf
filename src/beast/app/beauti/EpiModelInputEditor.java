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

package beast.app.beauti;

import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import epiinf.models.EpidemicModel;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.table.DefaultTableModel;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiModelInputEditor extends InputEditor.Base {

    EpidemicModel epidemicModel;
    String[] modelNames = {"Linear birth-death", "SIS", "SIR"};

    ComboBoxModel<String> emSelectorModel;
    DefaultTableModel psiSamplingRateModel, getPsiSamplingRateChangeTimesModel;

    public EpiModelInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return EpidemicModel.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr,
                     ExpandOption bExpandOption, boolean bAddButtons) {

        // Set up fields
        m_bAddButtons = bAddButtons;
        m_input = input;
        m_plugin = plugin;
        this.itemNr = itemNr;
        epidemicModel = (EpidemicModel) input.get();

        // Adds label to the left of input editor
        addInputLabel();

        // Create component models and fill them with data from input
        emSelectorModel = new DefaultComboBoxModel<>(modelNames);

        // Create and lay out GUI components
        JPanel panel = new JPanel(new GridBagLayout());
        panel.setBorder(new EtchedBorder());

        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(3, 3, 3, 3);
        c.weighty = 0.5;

        c.gridx = 0;
        c.gridy = 0;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Epidemic model to use:"), c);

        JComboBox<String> emSelector = new JComboBox<>(emSelectorModel);
        c.gridx = 1;
        c.gridy = 0;
        c.weightx = 1.0;
        c.anchor = GridBagConstraints.LINE_START;
        panel.add(emSelector, c);

        c.gridx = 0;
        c.gridy = 1;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Psi sampling rate:"), c);

        c.gridx = 0;
        c.gridy = 2;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JCheckBox("Psi sampling change times:"), c);

        c.gridx = 0;
        c.gridy = 3;
        c.weightx = 0.0;
        c.anchor = GridBagConstraints.LINE_END;
        panel.add(new JLabel("Rho sampling prob:"), c);

        add(panel);

        // Event handlers
    }

    class EpiParamModel {
        int nChanges;
    }
}
