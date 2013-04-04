package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:32
 * To change this template use File | Settings | File Templates.
 */
public class GaussianApodizator extends AbstractApodizator {

    public GaussianApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    public GaussianApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
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

        double time = i*processing.getDwellTime();
        return Math.exp(-2*Math.pow(time*lbGauss,2));
    }

}
