/*
 * Copyright (C) 2014 Tim Vaughan (tgvaughan@gmail.com)
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
import beast.util.TreeParser;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

/**
 * @author Tim Vaughan (tgvaughan@gmail.com)
 */
public class TreeParserFromFile extends TreeParser {

    public Input<String> newickFileNameInput = new Input<>("newickFileName",
        "Name of file containing tree in Newick format.");

    public TreeParserFromFile() { }

    @Override
    public void initAndValidate() throws Exception {
        if (newickFileNameInput.get() != null) {
            FileInputStream inputStream = new FileInputStream(newickFileNameInput.get());
            BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
            newickInput.setValue(reader.readLine(), this);
        }
        super.initAndValidate();
    }

    
    
}
