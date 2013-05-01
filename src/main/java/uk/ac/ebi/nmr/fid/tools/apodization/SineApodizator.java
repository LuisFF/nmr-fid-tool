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
 * Applies a Sine Bell window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 15:18
 *
 */
public class SineApodizator extends AbstractApodizator {

    protected SineApodizator(Spectrum spectrum) {
        super(spectrum);
    }
    //TODO fix this class
    public Spectrum calculate() {

        double factor;
        // ssbSine is defined as 180/ssb
        // TODO allow to redifine
        double offset = (180.0 - spectrum.getProc().getSsbSine())/180.0;
        // check the loop conditions in the original code...
        for (int i =0 ;
             i< (spectrum.getProc().getTdEffective()-1) &&
                     i >= ((spectrum.getProc().getTdEffective()-1)-spectrum.getRealChannelData().length*2) ;
             i+=2){
            factor=Math.sin(i/spectrum.getProc().getTdEffective()*Math.PI*offset);
            spectrum.setRealChannelData(spectrum.getProc().getTdEffective()-1 -i,
                    spectrum.getRealChannelData()[spectrum.getProc().getTdEffective()-1 -i]*factor);
            spectrum.setImaginaryChannelData(spectrum.getProc().getTdEffective()-2 -i,
                    spectrum.getImaginaryChannelData()[spectrum.getProc().getTdEffective()-2 -i]*factor);
        }
        return spectrum;
    }

    @Override
    protected double calculateFactor(int i, double lineBroadening) throws Exception {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected double calculateFactor(int i) throws Exception {
        double offset = (180.0 - spectrum.getProc().getSsbSine())/180.0;

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
