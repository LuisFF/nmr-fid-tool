package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 15:16
 * To change this template use File | Settings | File Templates.
 */
public class ZeroFillApodizator extends AbstractApodizator {
    public ZeroFillApodizator(double[] spectrum, Proc processing) {
        super(spectrum, processing);
    }

    public ZeroFillApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        super(spectrum, acquisitionMode, processing);
    }

    public double[] calculate(){
        return run((int) Math.floor(spectrum.length/5));
    }

    /**
     * This operation adds a determined set of zeros to the fid tail.
     * @param numberZeros
     */
    private double[] run(int numberZeros) {
        double[] data;
        //TODO Check if the sequencial and simultaneous stuff is right...
        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            data= new double[spectrum.length+numberZeros];
        }   else {
            data= new double[spectrum.length+numberZeros*2];
        }
        System.arraycopy(spectrum, 0, data, 0, spectrum.length);
        for (int i=spectrum.length; i < data.length; i++)
            data[i]=0;
        return data;
    }
}
