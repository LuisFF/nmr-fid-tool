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
 * Applies a Lorentz to Gauss window function to the fid.
 *
 * @author Luis F. de Figueiredo
 *
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:25
 *
 */

public class LorentzGausApodizatior extends AbstractApodizator {
    protected LorentzGausApodizatior(Spectrum spectrum) {
        super(spectrum);
    }


    /**
     * calculates the weigth for the Lorentz-Guassian apodization: W(i)= (lb*(i*dw)*pi(1-t/(2*tmax)))
     * where t=i*dw gives the time coordinate in (s), lb is the line broadening conversion from Lorents to Gaussian
     * in Hz, the gb is the GB-factor ??? from in the Proc file, and N*dw is the acquisition duration in s.
     * The cuteNMR implementation is W(i) = exp(-lb * pi * t (1-t/(gb*2*(N*dw)))
     *
     *  I think there is an error in the cuteNMR implementation because the exponential term has to be multiplied by
     *  -1 or or at least have one of the constants negative (gb or lb)... otherwise the weighting function has higher
     *  values in the region with more withe noise...
     *
     * This weighting function seems to be called Gaussian resolution enhancement (GRE) function in Pearson (1987)
     * J. Mag. Res. 74: 541-545. (or perhpas this is the actual Gaussian weigthing function)
     *
     * In Pearson 1987: W(i)=exp(lb*(i*dw)*pi(1-t/(2*tmax)))
     * where tmax = (2ln (2))/(pi*lb*ro^2) = GB / AQ. The ro is the relative linewidth [(lb_0+ lb_LB)/lb_0 -- lb_0 is
     * the linebroadening and lb_LB is the Lorentzian linebroadening], the AQ is the acquisition time N*dw, and in
     * Bruker spectrometers GB corresponds to the GB factor in the processing files.
     *
     * @param i
     * @param lbLorentzToGauss
     * @return
     */

    @Override
    public double calculateFactor (int i, double lbLorentzToGauss) throws Exception {
        //TODO check if the equation to calculate the factor is correct
        if(lbLorentzToGauss>=0)
            throw new Exception ("the Lorentz to Gaussian linebroadening factor has to be negative");

        if(spectrum.getProc().getGbFactor()<= 0 || spectrum.getProc().getGbFactor() > 1)
            throw new Exception ("the GB-factor has to be in the interval ]0,1[");
        /*
        cuteNMR implementation
         */
        double acquisitionTime = spectrum.getProc().getDwellTime()*spectrum.getProc().getTdEffective();
        double a = Math.PI*lbLorentzToGauss;
        double b = -a/(2.0*spectrum.getProc().getGbFactor()*acquisitionTime);
        double time=i*spectrum.getProc().getDwellTime();
        return Math.exp(-a*time -b*time*time);
        /*
        Pearson GRE function
         */
//        double a = Math.PI*lbLorentzToGauss;
//        double tmax = processing.getTdEffective()*processing.getDwellTime()/processing.getGbFactor();
//        return Math.exp(a*i*processing.getDwellTime()*
//                (1-tmax*i*processing.getDwellTime()/2));
    }

    /**
     * calculates the Lourentz-Gaussian weight with the default Lorentz to Gauss line broadening conversion.
     *
     * @param i
     * @return
     */
    @Override
    public double calculateFactor (int i) throws Exception {
        return calculateFactor(i, -1);
    }
}
