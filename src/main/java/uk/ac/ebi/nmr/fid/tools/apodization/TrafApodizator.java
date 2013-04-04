package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 12:29
 * To change this template use File | Settings | File Templates.
 */

public class TrafApodizator extends AbstractApodizator {

    protected TrafApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    protected TrafApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
    }

    @Override
    protected double calculateFactor(int i, double lbTraf) throws Exception {
        double acquisitionTime = processing.getDwellTime()*processing.getTdEffective();
        double time, e , eps;
        time=i*processing.getDwellTime();
        e = Math.exp(-time*lbTraf*Math.PI);
        eps = Math.exp((time-acquisitionTime) * lbTraf * Math.PI);
        return e/(e*e+eps*eps);
    }

    @Override
    protected double calculateFactor(int i) throws Exception {
        return calculateFactor(i, 0.2);
    }

}
