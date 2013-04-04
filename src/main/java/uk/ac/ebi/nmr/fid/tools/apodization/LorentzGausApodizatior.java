package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:25
 * To change this template use File | Settings | File Templates.
 */

public class LorentzGausApodizatior extends AbstractApodizator {
    protected LorentzGausApodizatior(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    protected LorentzGausApodizatior(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
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

        if(processing.getGbFactor()<= 0 || processing.getGbFactor() > 1)
            throw new Exception ("the GB-factor has to be in the interval ]0,1[");
        /*
        cuteNMR implementation
         */
        double acquisitionTime = processing.getDwellTime()*processing.getTdEffective();
        double a = Math.PI*lbLorentzToGauss;
        double b = -a/(2.0*processing.getGbFactor()*acquisitionTime);
        double time=i*processing.getDwellTime();
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
