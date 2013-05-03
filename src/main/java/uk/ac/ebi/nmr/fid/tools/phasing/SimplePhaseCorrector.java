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

package uk.ac.ebi.nmr.fid.tools.phasing;

import uk.ac.ebi.nmr.fid.Spectrum;

/**
 * Method to perform phase correction
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/05/2013
 * Time: 14:47
 *
 */
public class SimplePhaseCorrector implements PhaseCorrector {

    @Override
    public Spectrum phaseCorrection(Spectrum spectrum, double zeroPhase, double firstOrderPhase) {
        return phaseCorrection(spectrum, zeroPhase, firstOrderPhase,0);
    }

    @Override
    public Spectrum phaseCorrection(Spectrum spectrum, double zeroPhase, double firstOrderPhase, int pivot) {
        double [] realChannel = new double[spectrum.getAcqu().getAquiredPoints()/2];
        double [] imaginaryChannel = new double[spectrum.getAcqu().getAquiredPoints()/2];

        if(pivot > spectrum.getAcqu().getAquiredPoints()/2 || pivot < 0 ) {
            try {
                throw new Exception (" pivot is incorrectly set");
            } catch (Exception e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        for (int i = 0; i<spectrum.getAcqu().getAquiredPoints()/2 ; i++){
            double phaseAngle = 2 * Math.PI / 360 * ( i / (spectrum.getFid().length/2) * firstOrderPhase + zeroPhase );
            realChannel[i]=spectrum.getRealChannelData()[i]*Math.cos(phaseAngle)
                    - spectrum.getImaginaryChannelData()[i]*Math.sin(phaseAngle);
            imaginaryChannel[i]=spectrum.getRealChannelData()[i]*Math.sin(phaseAngle)
                    + spectrum.getImaginaryChannelData()[i]*Math.cos(phaseAngle);
        }

        spectrum.setRealChannelData(realChannel);
        spectrum.setImaginaryChannelData(imaginaryChannel);
        return spectrum;
    }
}