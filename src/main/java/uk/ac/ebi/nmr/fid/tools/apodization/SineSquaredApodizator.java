package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 15:22
 * To change this template use File | Settings | File Templates.
 */
public class SineSquaredApodizator extends AbstractApodizator {

    protected SineSquaredApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    protected SineSquaredApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
    }

    public double [] calculate() {
        double [] data = new double[processing.getTdEffective()];
        double factor;
        double offset = (180.0 - processing.getSsbSineSquared())/180.0;
        // check the loop conditions in the original code...
        for (int i =0 ;
             i< (processing.getTdEffective()-1) && i >= ((processing.getTdEffective()-1)-spectrum.length) ;
             i+=2){
            factor=Math.pow(Math.sin(i/processing.getTdEffective()*Math.PI*offset),2);
            data[processing.getTdEffective()-1 -i]=spectrum[processing.getTdEffective()-1 -i]*factor;
            data[processing.getTdEffective()-2 -i]=spectrum[processing.getTdEffective()-2 -i]*factor;
        }
        return data;
    }
}
