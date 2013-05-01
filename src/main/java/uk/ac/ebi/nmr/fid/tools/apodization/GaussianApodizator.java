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
 * Applies a gaussian window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:32
 *
 */
public class GaussianApodizator extends AbstractApodizator {

    public GaussianApodizator(Spectrum spectrum) {
        super(spectrum);
    }


    /**
     * calculates the weigth for the guassian apodization: W(i)=exp(-2(i*dw*lb)^2)
     * where i*dw gives the time coordinate in (s) and sigma is the line broadening with default value 0.1 Hz.
     * The cuteNMR implementation is de facto W(i)=exp(-2(i*dw*lb)^2))
     * In Vogt 2004: W(i)=exp(-lb*(i*dw)^2))
     * In Wikipedia: W(i)=exp(-1/2*(((n-(N-1)/2))/(sigma * (N-1)/2))^2))
     *
     * @throws Exception
     */
    @Override
    protected double calculateFactor(int i) throws Exception {
        spectrum.getProc().setLineBroadening(0.1);
        return calculateFactor(i,0.1);
    }

    /**
     * performs the guassian apodization with specified line broadning.
     * @param lbGauss
     * @throws Exception
     */
    @Override
    protected double calculateFactor(int i, double lbGauss) throws Exception {

        if(lbGauss ==0)
            throw new Exception ("line broadening cannot be zero");
        // original expression: (1.0/std::sqrt(2))/par->lbGauss * (1.0/std::sqrt(2))/par->lbGauss;
//        double sigmaFactor=0.5*Math.pow(1/lbGauss,2);

        double time = i*spectrum.getProc().getDwellTime();
        return Math.exp(-2*Math.pow(time*lbGauss,2));
    }

}
