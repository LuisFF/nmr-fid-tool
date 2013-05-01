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
 * Applies a TRAF window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:30
 *
 */
public class TrafsApodizator extends AbstractApodizator{

    public TrafsApodizator(Spectrum spectrum) {
        super(spectrum);
    }

    @Override
    protected double calculateFactor(int i, double lbTraf) throws Exception {
        double acquisitionTime = spectrum.getProc().getDwellTime()*spectrum.getProc().getTdEffective();
        double time, e , eps;
        time=i*spectrum.getProc().getDwellTime();
        e = Math.exp(-time*lbTraf*Math.PI);
        eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
        return (e*e*(e+eps))/(e*e*e+eps*eps*eps);

    }

    @Override
    protected double calculateFactor(int i) throws Exception {
        //this creates imutability issues
//        spectrum.getProc().setLineBroadening(0.03);
        return calculateFactor(i, 0.03);
    }

}
