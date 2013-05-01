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

package uk.ac.ebi.nmr.fid.tools;

import uk.ac.ebi.nmr.fid.Proc;

/**
 * Performs phase correction to the transformed spectrum
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 18/02/2013
 * Time: 05:29
 *
 */
public class PhaseCorrectionTool {

    double [] spectrum;
    Proc processing;
    double [] data;

    public PhaseCorrectionTool(double[] spectrum, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
        this.data = new double[processing.getTdEffective()/2+1];
    }

    public double [] zeroOrderPhasing (double angle){
        for(int i =0; i < processing.getTdEffective()-1 ; i+=2){
            data[i/2]=spectrum[i]*Math.cos(angle)-spectrum[i+1]*Math.sin(angle);
        }
        return data;
    }

    public double [] firstOrderPhasing (double angle){
        for(int i=0; i< processing.getTdEffective()-1; i+=2){
            data[i/2]=spectrum[i]*Math.cos(i/2*angle/(processing.getTdEffective()/2))
                    +spectrum[i+1]*Math.sin(i/2*angle/(processing.getTdEffective()/2));
        }
        return data;
    }




}
