package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 11:41
 * To change this template use File | Settings | File Templates.
 */
public class ExponentialApodizator extends AbstractApodizator {

    public ExponentialApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    public ExponentialApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
    }

    /**
     * calculates the weigth for the guassian apodization: W(i)= exp(-i*dw*lb*Pi)
     * where i*dw gives the time coordinate in (s) and the line broadening (as defined in the acqu) is in Hz.
     * The cuteNMR implementation is W(i) = exp(-(dw*i) * lb * pi)
     * In Vogt 2004: W(i)=exp(-lb*(i*dw)^2))
     *
     * @param i
     * @return
     */
    @Override
    protected double calculateFactor(int i) {
        return Math.exp(-i*processing.getDwellTime()*processing.getLineBroadening()*Math.PI);
    }

    /**
     * performs the exponential apodization with specified line broadening.
     *
     * @param lineBroadening
     * @throws Exception
     */
    @Override
    protected double calculateFactor(int i, double lineBroadening) {
        return Math.exp(-i*processing.getDwellTime()*lineBroadening*Math.PI);
    }


    
    
}
