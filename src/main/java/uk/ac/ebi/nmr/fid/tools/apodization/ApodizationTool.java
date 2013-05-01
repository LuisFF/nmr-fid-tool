/**
 * 
 * Copyright (C) 2010  Pascal Fricke
 * Copyright (C) 2000-2010  Kirk Marat, The University of Manitoba
 * 
 * 
 * This file is part of nmr-fid library. 
 * nmr-fid library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * nmr-fid library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with cuteNMR.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * This code is based upon the libraries released by Dr Kirk Marat from the University
 * of Manitoba and cuteNMR:
 *      ftp://davinci.chem.umanitoba.ca/pub/marat/SpinWorks/source_library/
 *      https://sourceforge.net/projects/cutenmr/
 */
package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * The ApodizationTool allows to apply various methods to improve the signal to noise value of an fid.
 * 
 * Copyright (C) 2010  Pascal Fricke
 * Copyright (C) 2000-2010  Kirk Marat, The University of Manitoba
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 14/01/2013
 * Time: 13:54
 */

@Deprecated
public class ApodizationTool {
    private static double[] spectrum;
    private Acqu.AcquisitionMode acquisitionMode= Acqu.AcquisitionMode.SEQUENTIAL;
    private Proc processing;

    /**
     *
     * @param spectrum
     * @param processing
     */
    public ApodizationTool(double[] spectrum, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
    }
    // TODO perhaps change this to accept an Acqu object, to make aparent that the acquisition mode has to be defined
    public ApodizationTool(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        this.spectrum=spectrum;
        this.processing=processing;
        this.acquisitionMode=acquisitionMode;
    }

//    public double[] getSpectrum() {
//        return spectrum;
//    }

    /**
     * performs the exponential apodization of the fid using the function F(x,i)= x.exp(-i*dw*lb)
     * where i*dw gives the time coordinate in (s) and the line broadening is in Hz.
     * The cuteNMR implementation is W(x,i) = x.exp(-(dw*i) * lb * pi)
     * @return
     */
//    public double [] exponential(){
//        double data[] = new double[processing.getTdEffective()];
//        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
//            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++)
//                data[i]*= Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
////                data[i] = spectrum[i]* Math.exp(-i*processing.getDwellTime()*processing.getLineBroadning());
//        } else { // simultaneous data
//            for (int i =0 ; i< processing.getTdEffective()-1 && i < spectrum.length; i+=2){
//                double factor = Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
////                double factor = Math.exp(-i*processing.getDwellTime()*processing.getLineBroadning());
//                data[i]= spectrum[i]* factor;
//                data[i+1]=spectrum[i+1]* factor;
////                spectrum[i+1]*= Math.exp(-i*processing.getDwellTime()*Math.PI*processing.getLineBroadning());
//            }
//        }
//        return data;
//    }













    /**
     * performs the guassian apodization according to Bramer 2001: F(x,i)=x*exp(-(i*dw*lb)^2)
     * where i*dw gives the time coordinate in (s) and the line broadening (lbGauss) in Hz.
     * This also correspondes to one of the alternative ways of defining the standard gaussian distribution,
     * according to wikipedia.
     * @param lbGauss
     * @throws Exception
     */
    public double [] gaussianBramer2001(double lbGauss) throws Exception {
        double [] data = new double[processing.getTdEffective()];
        double time;
        if(lbGauss ==0)
            throw new Exception ("line broadening cannot be zero");
//        double sigmaFactor=0.5*Math.pow(1/lbGauss,2);

        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                time=i*processing.getDwellTime();
                data[i]=spectrum[i]*Math.exp(-Math.pow(time*lbGauss,2));
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                time=i*processing.getDwellTime();
                double factor=Math.exp(-Math.pow(time*lbGauss,2));
                data[i]=spectrum[i]*factor;
                data[i+1]=spectrum[i+1]*factor;
            }
        }
        return data;
    }

    public void firstPoint(){
        firstPoint(1);
    }
    public void firstPoint(double fpCorrection){
        spectrum[0]*=fpCorrection;
        if(!acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            spectrum[1]*=fpCorrection;
        }
    }


}
