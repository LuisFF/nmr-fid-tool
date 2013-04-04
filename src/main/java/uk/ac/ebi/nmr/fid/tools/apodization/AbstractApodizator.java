package uk.ac.ebi.nmr.fid.tools.apodization;

import uk.ac.ebi.nmr.fid.Acqu;
import uk.ac.ebi.nmr.fid.Proc;

/**
 * Created with IntelliJ IDEA.
 * User: ldpf
 * Date: 03/04/2013
 * Time: 11:42
 * To change this template use File | Settings | File Templates.
 */
abstract public class AbstractApodizator implements Apodizator {

    double[] spectrum;
    Acqu.AcquisitionMode acquisitionMode= Acqu.AcquisitionMode.SEQUENTIAL;
    Proc processing;

    public AbstractApodizator(double [] spectrum, Proc processing) {
        this.spectrum = spectrum;
        this.processing = processing;
    }

    // TODO perhaps change this to accept an Acqu object, to make apparent that the acquisition mode has to be defined
    public AbstractApodizator(double[] spectrum, Acqu.AcquisitionMode acquisitionMode, Proc processing) {
        this.spectrum = spectrum;
        this.acquisitionMode = acquisitionMode;
        this.processing = processing;
    }

    public double[] calculate() throws Exception {

        double [] data = new double[processing.getTdEffective()];
        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                data[i]=spectrum[i]*calculateFactor(i);
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                data[i]=spectrum[i]*calculateFactor(i);
                data[i+1]=spectrum[i+1]*calculateFactor(i);
            }
        }
        return data;
    }


    public double[] calculate(double lineBroadning) throws Exception {
        double [] data = new double[processing.getTdEffective()];
        if (acquisitionMode.equals(Acqu.AcquisitionMode.SEQUENTIAL)) {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i++){
                data[i]=spectrum[i]*calculateFactor(i,lineBroadning);
            }
        } else {
            for (int i =0 ; i< processing.getTdEffective() && i < spectrum.length; i+=2){
                data[i]=spectrum[i]*calculateFactor(i,lineBroadning);
                data[i+1]=spectrum[i+1]*calculateFactor(i, lineBroadning);
            }
        }
        return data;
    }


    protected abstract double calculateFactor(int i, double lineBroadening) throws Exception;

    protected abstract double calculateFactor(int i) throws Exception;
}
