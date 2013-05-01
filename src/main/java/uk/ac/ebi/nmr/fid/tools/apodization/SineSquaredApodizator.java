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

package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Spectrum;

/**
 *
 * Applies a Sine Squared window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 15:22
 */
public class SineSquaredApodizator extends AbstractApodizator {

    protected SineSquaredApodizator(Spectrum spectrum) {
        super(spectrum);
    }

    public Spectrum calculate() {

        double factor;

        // check the loop conditions in the original code...
        for (int i =0 ;
             i< (spectrum.getProc().getTdEffective()-1) &&
                     i >= ((spectrum.getProc().getTdEffective()-1)-spectrum.getRealChannelData().length*2) ;
             i+=2){
            try {
                spectrum.setImaginaryChannelData(spectrum.getProc().getTdEffective()-2 -i,
                        spectrum.getRealChannelData()[spectrum.getProc().getTdEffective()-2 -i]*calculateFactor(i));
                spectrum.setRealChannelData(spectrum.getProc().getTdEffective()-1 -i,
                        spectrum.getRealChannelData()[spectrum.getProc().getTdEffective()-1 -i]*calculateFactor(i));
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }
        return spectrum;
    }

    @Override
    protected double calculateFactor(int i, double lineBroadening) throws Exception {
        return calculateFactor(i);
    }

    @Override
    protected double calculateFactor(int i) throws Exception {
        double offset = (180.0 - spectrum.getProc().getSsbSineSquared())/180.0;
        return Math.pow(Math.sin(i/spectrum.getProc().getTdEffective()*Math.PI*offset),2);
    }
}
