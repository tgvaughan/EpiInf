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

/**
 * Basic tuple.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 * @param <O>
 * @param <T>
 * @param <X>
 */
public class Triple<O,T,X> {

    public O one;
    public T two;
    public X three;

    public Triple(O one, T two, X three) {
        this.one = one;
        this.two = two;
        this.three = three;
    }
}
