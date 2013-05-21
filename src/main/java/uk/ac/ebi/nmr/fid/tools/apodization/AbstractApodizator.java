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

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Spectrum;

/**
 * Abstract class that applies a window function to the fid in order to reduce noise or enhance signal.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 11:42
 *
 */
abstract public class AbstractApodizator implements Apodizator {

    Spectrum spectrum;


    public AbstractApodizator(Spectrum spectrum) {
        this.spectrum = spectrum;
    }

    public Spectrum calculate() throws Exception {
        // perhaps clone the spectrum otherwise one can have more clashes with the acqu and proc objects
        // because of imutability issues
        double [] fid = new double[spectrum.getFid().length];
        if (spectrum.getAcqu().getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< spectrum.getProc().getTdEffective() && i < spectrum.getFid().length; i++){
                fid[i]=spectrum.getRealChannelData()[i] * calculateFactor(i);
            }
        } else {
            for (int i =0 ; i < spectrum.getFid().length/2; i++){
                // set real values, i.e. even numbers
                fid[i*2] = spectrum.getFid()[i * 2] * calculateFactor(i);
                // set imaginary values, i.e. odd numbers
                fid[i * 2 + 1]= spectrum.getFid()[i * 2 + 1] * calculateFactor(i);
            }
        }
        return new Spectrum(fid,spectrum.getAcqu(),spectrum.getProc());
    }

    public Spectrum calculate(double lineBroadning) throws Exception {
        double [] fid = new double[spectrum.getFid().length];
        // this creates issues with immutability, need to clone the proc object
//        spectrum.getProc().setLineBroadening(lineBroadning);
        if (spectrum.getAcqu().getAcquisitionMode().equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< spectrum.getProc().getTdEffective() && i < spectrum.getFid().length; i++){
                fid[i]=spectrum.getFid()[i] * calculateFactor(i, lineBroadning);
            }
        } else {
            for (int i =0 ; i < spectrum.getFid().length/2; i++){
                // set real values, i.e. even numbers
                fid[i * 2]= spectrum.getFid()[i * 2] * calculateFactor(i, lineBroadning);
                // set imaginary values, i.e. odd numbers
                fid[i * 2 + 1]= spectrum.getFid()[i * 2 + 1] * calculateFactor(i, lineBroadning);
            }
        }
        return new Spectrum(fid,spectrum.getAcqu(),spectrum.getProc());
    }


    protected abstract double calculateFactor(int i, double lineBroadening) throws Exception;

    protected abstract double calculateFactor(int i) throws Exception;
}
