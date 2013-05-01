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
 * Applies an exponential window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 11:41
 *
 */
public class ExponentialApodizator extends AbstractApodizator {

    public ExponentialApodizator(Spectrum spectrum) {
        super(spectrum);
    }

    /**
     * calculates the weigth for the exponential apodization: W(i)= exp(-i*dw*lb*Pi)
     * where i*dw gives the time coordinate in (s) and the line broadening (as defined in the acqu) is in Hz.
     * The cuteNMR implementation is W(i) = exp(-(dw*i) * lb * pi)
     * In Vogt 2004: W(i)=exp(-lb*(i*dw)^2))
     *
     * @param i
     * @return
     */
    @Override
    protected double calculateFactor(int i) {
        return Math.exp(-i*spectrum.getProc().getDwellTime()*spectrum.getProc().getLineBroadening()*Math.PI);
    }

    /**
     * performs the exponential apodization with specified line broadening.
     *
     * @param lineBroadening
     * @throws Exception
     */
    @Override
    protected double calculateFactor(int i, double lineBroadening) {
        return Math.exp(-i*spectrum.getProc().getDwellTime()*lineBroadening*Math.PI);
    }


    
    
}
