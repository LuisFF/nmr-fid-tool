/*
 * Copyright (c) 2013. EMBL, European Bioinformatics Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package uk.ac.ebi.nmr.fid;

import org.apache.commons.lang3.ArrayUtils;

import java.util.List;

/**
 * Data structure for the fid
 *
 * @author  Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 14/01/2013
 * Time: 14:00
 *
 */

@Deprecated
public class Fid {
    private double[] data;
    private double[] real;
    private double[] imaginary;

    public Fid(List<Integer> fid) {
        this(ArrayUtils.toPrimitive(fid.toArray(new Integer[fid.size()])));
    }

    public Fid(double [] fid) {
        data = fid;
        splitData();

    }

    private void splitData() {
        real=new double[data.length/2];
        imaginary=new double [data.length/2];
        for(int i=0; i < data.length; i+=2){
            real[i/2]=data[i];// real are in even positions
            imaginary[i/2]=data[i+1];// imaginary are in odd positions
        }
    }

    public Fid(int [] fid) {
        data = new double[fid.length];
        for(int i = 0 ; i < fid.length; i++){
            data[i]= (double) fid[i];
        }
        splitData();
    }

    public Fid(float [] fid) {
        data = new double[fid.length];
        for(int i = 0 ; i < fid.length; i++){
            data[i]= (double) fid[i];
        }
        splitData();
    }

    public double[] getData() {
        return data;
    }

    public double[] getReal() {
        return real;
    }

    public double[] getImaginary() {
        return imaginary;
    }
}
