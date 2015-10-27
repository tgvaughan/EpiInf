/*
 * Copyright (C) 2014 Tim Vaughan <tgvaughan@gmail.com>
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

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import com.google.common.collect.Lists;
import epiinf.TreeEvent;
import epiinf.TreeEventList;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.List;

/**
 * Some utility methods that don't fit anywhere else.
 *
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public class EpiInfUtilityMethods {

    /**
     * Translate event list into file readable by expoTree's calc_likelihood.
     * 
     * @param tree TreeEventList for a tree.
     * @param origin Time of origin.
     * @param ps PrintStream  where output is sent.
     */
    public static void writeExpoTreeFile(Tree tree, double origin, PrintStream ps) {

        List<Node> nodeList = Lists.newArrayList(tree.getNodesAsArray());
        nodeList.sort((Node o1, Node o2) -> {
            if (o1.getHeight()<o2.getHeight())
                return -1;
            if (o1.getHeight()>o2.getHeight())
                return 1;
            return 0;
        });

        for (Node node : nodeList) {
            if (node.isLeaf())  {
                if (node.getHeight()>0.0)
                    ps.println(node.getHeight() + " 0");
            } else {
                ps.println(node.getHeight() + " 1");
            }
        }

        ps.println(origin + " 99");
    }
    
}
