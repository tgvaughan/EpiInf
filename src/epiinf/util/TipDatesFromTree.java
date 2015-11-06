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

package epiinf.util;

import beast.core.Input;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class TipDatesFromTree extends TraitSet {
    public Input<Tree> treeInput = new Input<>("tree", "Tree from which to " +
            "extract tip dates.", Input.Validate.REQUIRED);

    public TipDatesFromTree() {
        traitsInput.setRule(Input.Validate.OPTIONAL);
        traitNameInput.setRule(Input.Validate.OPTIONAL);
        taxaInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() throws Exception {

        traitNameInput.setValue("date-backward", this);

        StringBuilder valueBuilder = new StringBuilder();

        boolean isFirst = true;
        for (Node leaf : treeInput.get().getExternalNodes()) {
            if (isFirst)
                isFirst = false;
            else
                valueBuilder.append(",");
            valueBuilder.append(leaf.getID() + "=" + leaf.getHeight());
        }
        traitsInput.setValue(valueBuilder.toString(), this);

        // Build taxon set if not available
        if (taxaInput.get() == null) {
            List<Taxon> taxonList = new ArrayList<>();
            for (Node leaf : treeInput.get().getExternalNodes())
                taxonList.add(new Taxon(leaf.getID()));
            taxaInput.setValue(new TaxonSet(taxonList), this);
        }

        super.initAndValidate();
    }
}
