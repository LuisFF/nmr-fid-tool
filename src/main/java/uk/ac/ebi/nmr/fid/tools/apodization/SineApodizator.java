package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 15:18
 * To change this template use File | Settings | File Templates.
 */
public class SineApodizator extends AbstractApodizator {

    protected SineApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    protected SineApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
    }

    public double [] calculate() {
        double [] data = new double[processing.getTdEffective()];
        double factor;
        // ssbSine is defined as 180/ssb
        //
        double offset = (180.0 - processing.getSsbSine())/180.0;
        // check the loop conditions in the original code...
        for (int i =0 ;
             i< (processing.getTdEffective()-1) && i >= ((processing.getTdEffective()-1)-spectrum.length) ;
             i+=2){
            factor=Math.sin(i/processing.getTdEffective()*Math.PI*offset);
            data[processing.getTdEffective()-1 -i]=spectrum[processing.getTdEffective()-1 -i]*factor;
            data[processing.getTdEffective()-2 -i]=spectrum[processing.getTdEffective()-2 -i]*factor;
        }
        return data;
    }

    @Override
    protected double calculateFactor(int i, double lineBroadening) throws Exception {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected double calculateFactor(int i) throws Exception {
        double offset = (180.0 - processing.getSsbSine())/180.0;

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
